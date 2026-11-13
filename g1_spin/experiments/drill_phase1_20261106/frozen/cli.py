"""Headless demo commands; see docs/algorithm.md for evidence and limitations."""
import argparse
import csv
import json
from pathlib import Path
import time

import cv2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import yaml

from .events import save_npz
from .io import source, windows
from .synthetic import generate
from .baseline import estimate_periodicity
from .spherical_cmax import lift, warp, iwe, estimate


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False))


def figure(path):
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()


def inspect_file(path, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    meta, batches = source(path)
    total = on = duplicate = max_gap = 0
    first = last = None
    histogram = {}
    for b in batches:
        b.validate(meta.get('width'), meta.get('height'))
        if not len(b):
            continue
        if last is not None:
            if int(b.t_us[0]) < last:
                raise ValueError('Timestamp regression across batches')
            duplicate += int(int(b.t_us[0]) == last)
            max_gap = max(max_gap, int(b.t_us[0])-last)
        quality = b.validate()
        duplicate += quality['duplicate_timestamps']
        max_gap = max(max_gap, quality['max_gap_us'])
        first = int(b.t_us[0]) if first is None else first
        last = int(b.t_us[-1])
        total += len(b)
        on += int(b.polarity.sum())
        bins, counts = np.unique((b.t_us-first)//1000, return_counts=True)
        for k, count in zip(bins, counts):
            histogram[int(k)] = histogram.get(int(k), 0) + int(count)
    duration = (last-first)*1e-6 if total else 0
    result = dict(metadata=meta, events=total, start_us=first, end_us=last,
                  duration_s=duration, mean_events_s=total/duration if duration else None,
                  on_fraction=on/total if total else None, off_fraction=(total-on)/total if total else None,
                  timestamp_regressions=0, duplicate_timestamps=duplicate, max_gap_us=max_gap)
    write_json(output/'inspection.json', result)
    with (output/'event_rate.csv').open('w') as f:
        writer = csv.writer(f); writer.writerow(['relative_ms', 'events'])
        writer.writerows((k, histogram.get(k, 0)) for k in range(max(histogram, default=-1)+1))
    plot_bins = sorted(set(histogram) | {k-1 for k in histogram if k > 0} | {k+1 for k in histogram})
    plt.figure(figsize=(9, 3)); plt.plot(plot_bins, [histogram.get(k, 0) for k in plot_bins], lw=.8)
    plt.xlabel('Time (ms)'); plt.ylabel('Events / ms'); figure(output/'event_rate.png')
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def preview(path, output, window_ms=2, start_us=None, duration_ms=100, circle=None):
    output = Path(output); output.parent.mkdir(parents=True, exist_ok=True)
    meta, batches = source(path)
    width, height = meta.get('width') or 256, meta.get('height') or 256
    if meta.get('format') == 'g1_spin EventBatch NPZ':
        width, height = max(256, width), max(256, height)
    # MP4 encoders require even image sizes.
    width += width % 2; height += height % 2
    writer = None
    if output.suffix.lower() == '.mp4':
        writer = cv2.VideoWriter(str(output), cv2.VideoWriter_fourcc(*'mp4v'), 50, (width, height))
        if not writer.isOpened():
            raise RuntimeError('MP4 encoder failed')
    counts = np.zeros((height, width), np.int64)
    surface = np.full((height, width), -1., float)
    frames = 0
    origin_us = None
    try:
        for start, end, b in windows(batches, int(window_ms*1000), start_us,
                                     int(duration_ms*1000) if duration_ms else None):
            b.validate(width, height)
            x, y = np.floor(b.x).astype(int), np.floor(b.y).astype(int)
            frame = np.full((height, width, 3), 15, np.uint8)
            frame[y[b.polarity == 1], x[b.polarity == 1]] = [80, 200, 255]
            frame[y[b.polarity == 0], x[b.polarity == 0]] = [240, 130, 80]
            np.add.at(counts, (y, x), 1)
            origin_us = start if origin_us is None else origin_us
            np.maximum.at(surface, (y, x), (b.t_us-origin_us)*1e-6)
            if circle:
                cx, cy, radius, ratio = circle
                cv2.circle(frame, (round(cx), round(cy)), round(radius), (180, 180, 180), 1)
                cv2.circle(frame, (round(cx), round(cy)), round(radius*ratio), (80, 220, 100), 1)
            cv2.putText(frame, f'{start/1000:.1f} ms | {len(b)} events', (8, 18), 0, .4, (230, 230, 230), 1)
            if frames == 0:
                cv2.imwrite(str(output.with_suffix('.png')), frame)
            if writer:
                writer.write(frame)
            frames += 1
    finally:
        if writer:
            writer.release()
    if not frames:
        if writer:
            output.unlink(missing_ok=True)
        raise ValueError('No preview frames: select a time range containing events or set an explicit start and duration')
    plt.figure(figsize=(9, 4))
    plt.subplot(121); plt.imshow(np.log1p(counts), cmap='magma'); plt.title('Log event count')
    plt.subplot(122); plt.imshow(np.ma.masked_less(surface, 0), cmap='viridis'); plt.title('Latest event since preview start (s)')
    figure(output.with_name(output.stem+'_surfaces.png'))
    write_json(output.with_suffix('.json'), dict(frames=frames, fps=50, window_ms=window_ms,
               playback_speed=window_ms/20, timestamp_unit='us', source=str(path)))
    return frames


def analyze(path, output, config, truth=None):
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    ball, settings = config['ball'], config['analysis']
    cx, cy, radius = [ball.get(k) for k in ('center_x', 'center_y', 'radius_px')]
    ratio = ball.get('inner_radius_ratio', .9)
    if any(v is None for v in (cx, cy, radius)) or not np.isfinite([cx,cy,radius,ratio]).all() or radius <= 0 or not 0 < ratio < 1:
        raise ValueError('Set measured ball center_x, center_y, positive radius_px and 0<inner_radius_ratio<1')
    meta, batches = source(path)
    records = []; diagnostics = None
    window_us = round(settings.get('window_ms', 10)*1000)
    step_us = round(settings.get('step_ms', settings.get('window_ms',10))*1000)
    duration_ms = settings.get('duration_ms',100)
    started = time.perf_counter()
    for start, end, raw in windows(batches, window_us, settings.get('start_us'),
                                 round(duration_ms*1000) if duration_ms else None, step_us):
        raw.validate(meta.get('width'), meta.get('height'))
        mask = (raw.x-cx)**2 + (raw.y-cy)**2 < (radius*ratio)**2
        b = raw.select(mask)
        points = lift(b.x, b.y, cx, cy, radius)
        dt = b.seconds()
        before = time.perf_counter()
        result = estimate(points, dt, max_rpm=settings.get('max_rpm',3500))
        result.update(start_us=start, end_us=end, raw_events=len(raw), effective_events=len(b),
                      runtime_s=time.perf_counter()-before)
        result['quality']['outside_roi_fraction'] = 1-len(b)/len(raw) if len(raw) else None
        result['quality']['feature_visibility'] = 'not_independently_identified'
        if records and result['valid'] and records[-1]['valid']:
            previous = np.array(records[-1]['omega_rad_s'])
            if np.linalg.norm(np.array(result['omega_rad_s'])-previous) > .5*max(np.linalg.norm(previous),1):
                result['quality']['temporal_candidate_omega_rad_s'] = result['omega_rad_s']
                result.update(valid=False, failure_reason='adjacent_window_inconsistency',
                              omega_rad_s=None, axis_camera=None, spin_rpm=None, spin_rps=None, spin_rad_s=None)
        if diagnostics is None or (result['valid'] and not diagnostics[2]['valid']):
            diagnostics = (points, dt, result)
        records.append(result)
        print(f'window {start}:{end} events={len(b)} valid={result["valid"]} rpm={result["spin_rpm"]}', flush=True)
    good = [r for r in records if r['valid']]
    rpm = np.array([r['spin_rpm'] for r in good])
    summary = dict(source=str(path), synthetic_truth=truth, window_count=len(records),
                   valid_fraction=len(good)/len(records) if records else 0,
                   median_rpm=float(np.median(rpm)) if len(rpm) else None,
                   rpm_std=float(np.std(rpm)) if len(rpm) else None,
                   runtime_s=time.perf_counter()-started,
                   caveat='Research demo; idealized synthetic accuracy is not physical G1 accuracy.', windows=records)
    if truth and good:
        target = np.array(truth['axis'], float)*np.sign(truth['rpm']); target /= np.linalg.norm(target)
        summary['rpm_absolute_error'] = abs(summary['median_rpm']-abs(truth['rpm']))
        summary['rpm_relative_error'] = summary['rpm_absolute_error']/abs(truth['rpm'])
        summary['median_axis_error_deg'] = float(np.median([np.degrees(np.arccos(np.clip(np.dot(r['axis_camera'],target),-1,1))) for r in good]))
    write_json(output/'result.json', summary)
    columns = ['start_us','end_us','raw_events','effective_events','valid','spin_rpm','spin_rps','spin_rad_s','omega_rad_s','axis_camera','failure_reason','runtime_s']
    with (output/'windows.csv').open('w') as f:
        writer = csv.DictWriter(f,fieldnames=columns,extrasaction='ignore'); writer.writeheader(); writer.writerows(records)
    plt.figure(figsize=(9,3))
    plt.plot([(r['start_us']-records[0]['start_us'])/1000 for r in records], [r['spin_rpm'] if r['valid'] else np.nan for r in records], 'o-')
    if truth: plt.axhline(abs(truth['rpm']), color='gray', ls='--', label='Synthetic truth'); plt.legend()
    plt.xlabel('Window start (ms)'); plt.ylabel('RPM (valid windows only)'); figure(output/'rpm_curve.png')
    if diagnostics:
        p, dt, result = diagnostics
        omega = result['omega_rad_s'] if result['valid'] else [0,0,0]
        center_dt = dt-(dt[0]+dt[-1])/2 if len(dt) else dt
        images = [iwe(p), iwe(warp(p, center_dt, omega))]
        plt.figure(figsize=(10,4)); vmax=max(float(im.max()) for im in images)
        for i, im in enumerate(images):
            plt.subplot(1,2,i+1); plt.imshow(im,origin='lower',cmap='magma',vmin=0,vmax=vmax)
            plt.title(['Uncompensated spherical IWE','Compensated IWE' if result['valid'] else 'Invalid: no compensation claimed'][i]); plt.xlabel('Longitude'); plt.ylabel('Latitude')
        figure(output/'spherical_before_after.png')
        plt.figure(figsize=(10,3))
        plt.subplot(121)
        for h in result['quality'].get('objective_histories',[]): plt.plot(h,alpha=.5)
        plt.xlabel('Objective evaluation'); plt.ylabel('Normalized IWE variance')
        plt.subplot(122); scan=result['quality'].get('candidate_scan',[])
        plt.plot([s['rpm'] for s in scan],[s['contrast'] for s in scan]); plt.xlabel('Grid speed (RPM)'); plt.ylabel('Best axis contrast')
        figure(output/'optimization.png')
    preview(path, output/'event_preview.mp4',2,settings.get('start_us'),duration_ms,(cx,cy,radius,ratio))
    failures = {reason:sum(r['failure_reason']==reason for r in records) for reason in set(r['failure_reason'] for r in records) if reason}
    (output/'report.md').write_text('# 固定球三维旋转 Demo\n\n'
        + ('数据：理想球面特征轨迹合成事件，非物理 DVS 仿真。\n\n' if truth else '数据：输入文件；是否为合格旋转实验须人工确认。\n\n')
        + f'有效窗口：{len(good)}/{len(records)}；RPM 中位数：{summary["median_rpm"]}；标准差：{summary["rpm_std"]}。\n\n'
        + f'失败原因：{failures}。计算耗时：{summary["runtime_s"]:.2f} s。\n\n'
        + f'合成误差：RPM 相对误差 {summary.get("rpm_relative_error")}；有向轴误差 {summary.get("median_axis_error_deg")}°。\n\n'
        + '完整向量与质量诊断见 result.json，逐窗结果见 windows.csv。优化收敛不证明物理正确；多初值一致性和阈值是工程启发式。真实标记可观测性、热点、背景干扰仍需人工检查。\n')
    return summary


def benchmark(output):
    output=Path(output); output.mkdir(parents=True,exist_ok=True)
    cases=[('500rpm',500,(.3,.7,.6),0),('1000rpm',1000,(.3,.7,.6),0),
           ('2000rpm',2000,(-.5,.2,.8),0),('3000rpm',3000,(0,0,1),0),
           ('reverse',-1000,(.3,.7,.6),0),('background_noise',1000,(.3,.7,.6),.1),('empty',0,(0,0,1),0)]
    summaries=[]
    for name,rpm,axis,noise in cases:
        folder=output/name; folder.mkdir(exist_ok=True)
        b=generate(rpm=rpm,axis=axis,duration_s=.03,noise_fraction=noise)
        path=folder/'synthetic.npz'; save_npz(path,b)
        config=dict(ball=dict(center_x=128,center_y=128,radius_px=100,inner_radius_ratio=.9),
                    analysis=dict(window_ms=10,step_ms=10,start_us=0,duration_ms=30,max_rpm=3500))
        result=analyze(path,folder,config,dict(rpm=rpm,axis=axis,model='idealized_surface_landmarks'))
        result.pop('windows'); result['case']=name
        # Period baseline needs >=3 cycles: evaluate independent longer sample.
        long_sample = generate(rpm=rpm,axis=axis,duration_s=.6,noise_fraction=noise)
        tick = time.perf_counter()
        result['period_baseline']=estimate_periodicity(long_sample)
        result['period_baseline']['runtime_s'] = time.perf_counter()-tick
        result['period_baseline']['duration_s'] = .6
        try:
            from .optical_flow import estimate_flow
            tick = time.perf_counter()
            result['flow_baseline']=estimate_flow(b,128,128,100)
            result['flow_baseline']['runtime_s'] = time.perf_counter()-tick
            result['flow_baseline']['duration_s'] = .03
        except ImportError:
            result['flow_baseline']={'valid':False,'failure_reason':'not_implemented'}
        summaries.append(result)
    write_json(output/'benchmark.json',summaries)
    plt.figure(figsize=(10,4)); x=np.arange(len(summaries))
    plt.bar(x,[s['median_rpm'] or 0 for s in summaries],color='#289cae',label='Estimated (0 = invalid)')
    plt.scatter(x,[abs(c[1]) for c in cases],color='#ec8a32',label='Synthetic truth',zorder=3)
    plt.xticks(x,[s['case'] for s in summaries],rotation=20); plt.ylabel('RPM'); plt.legend(); figure(output/'benchmark.png')
    lines=['# 合成 Demo 实测结果','', '理想球面特征轨迹；没有向估计器传入真值或点对应。不是实际 G1 精度。','', '| Case | 有效窗口比例 | RPM | 相对误差 | 轴误差° |','|---|---:|---:|---:|---:|']
    for s in summaries:
        lines.append(f'| {s["case"]} | {s["valid_fraction"]:.2f} | {s["median_rpm"]} | {s.get("rpm_relative_error")} | {s.get("median_axis_error_deg")} |')
    lines += ['', '周期法只输出候选频率，无三维轴；10 ms 通常不足以覆盖多个周期，使用独立 0.6 s 片段。光流对照如有效仅为工程法向约束版本。细节见 benchmark.json。']
    (output/'report.md').write_text('\n'.join(lines))
    from .demo_report import write_demo
    write_demo(output, summaries)
    return summaries


def main():
    parser=argparse.ArgumentParser(description=__doc__); subs=parser.add_subparsers(dest='command',required=True)
    p=subs.add_parser('inspect'); p.add_argument('path'); p.add_argument('--output',default='outputs/inspect')
    p=subs.add_parser('preview'); p.add_argument('path'); p.add_argument('--output',default='outputs/preview.mp4'); p.add_argument('--window-ms',type=float,default=2); p.add_argument('--start-us',type=int); p.add_argument('--duration-ms',type=float,default=100)
    p=subs.add_parser('analyze'); p.add_argument('path'); p.add_argument('--config',default='configs/phase1.yaml'); p.add_argument('--output',default='outputs/analysis')
    p=subs.add_parser('synthetic-benchmark'); p.add_argument('--output',default='outputs/synthetic')
    args=parser.parse_args()
    if args.command=='inspect': inspect_file(args.path,args.output)
    elif args.command=='preview': preview(args.path,args.output,args.window_ms,args.start_us,args.duration_ms)
    elif args.command=='analyze': analyze(args.path,args.output,yaml.safe_load(Path(args.config).read_text()))
    else: benchmark(args.output)


if __name__=='__main__':
    main()
