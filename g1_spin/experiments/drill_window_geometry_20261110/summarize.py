"""Descriptive reporting; no tuning or unconfirmed cross-recording RPM scoring."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
rows=[json.loads(s) for s in (OUT/'results.jsonl').read_text().splitlines()]
assert len(rows)==126
methods=['fixed','linear','median3']
reference=json.loads((OUT/'reference.json').read_text())
nominal=reference['nominal_same_setting_rpm']
low,high=reference['same_setting_readout_envelope_rpm']
groups=[]
for run,condition in [(1,'rotation'),(4,'rotation'),(5,'translation'),(5,'static')]:
    for ms in [10,20,40]:
        for method in methods:
            group=[r for r in rows if (r['run'],r['condition'],r['window_ms'],r['method'])==
                   (run,condition,ms,method)]
            good=[r['rpm'] for r in group if r['valid']]
            groups.append(dict(run=run,condition=condition,window_ms=ms,method=method,
                planned=len(group),valid=len(good),median=float(np.median(good)) if good else None,
                minimum=min(good) if good else None,maximum=max(good) if good else None,
                approximate_reference_abs_deviation_mean_rpm=float(np.mean(np.abs(np.asarray(good)-nominal))) if good and condition=='rotation' else None,
                within_reference_envelope=sum(low<=v<=high for v in good) if condition=='rotation' else None))
(OUT/'by_run.json').write_text(json.dumps(groups,indent=2))

# Independently recover reference periods from saved marker coordinates and original PTS.
reference_checks=[]
for n in [1,2]:
    path=ROOT/f'experiments/wechat_spin_20261105/originals/video{n}_marker_A.csv'
    records=list(csv.DictReader(path.open()))
    crossings=[]
    for cycle in range(4):
        points=[r for r in records if int(r['cycle'])==cycle]
        found=[]
        for a,b in zip(points,points[1:]):
            x1,x2=float(a['marker_dx_over_ball_radius']),float(b['marker_dx_over_ball_radius'])
            if x1*x2<0:
                t1,t2=float(a['source_time_s']),float(b['source_time_s'])
                found.append((x2*t1-x1*t2)/(x2-x1))
        assert len(found)==1
        crossings.extend(found)
    turns=np.arange(4); times=np.asarray(crossings)
    period=np.sum((turns-turns.mean())*(times-times.mean()))/np.sum((turns-turns.mean())**2)
    saved=json.loads((OUT/'reference.json').read_text())['videos'][n-1]['marker_A_crossing_rpm']
    assert abs(60/period-saved)<1e-8
    reference_checks.append(dict(video=n,rpm=60/period,verified_from_marker_csv=True))
(OUT/'reference_verification.json').write_text(json.dumps(reference_checks,indent=2))

fig,axes=plt.subplots(2,2,figsize=(12,8))
for ax,(run,condition) in zip(axes.flat,[(1,'rotation'),(4,'rotation'),(5,'translation'),(5,'static')]):
    for method,offset,color in zip(methods,[-1.5,0,1.5],['#2878b5','#db7b20','#38965c']):
        group=[r for r in rows if r['run']==run and r['condition']==condition and r['method']==method and r['valid']]
        ax.scatter([r['window_ms']+offset for r in group],[r['rpm'] for r in group],label=method,color=color,s=35,alpha=.7)
    ax.set(title=f'run{run} {condition}: accepted candidates only',xlabel='Window (ms)',ylabel='RPM',xticks=[10,20,40])
    ax.grid(alpha=.2)
    if condition=='rotation':
        ax.axhspan(low,high,color='gray',alpha=.12)
        ax.axhline(nominal,color='gray',ls='--',lw=1,label='same-setting reference')
    ax.legend()
    if not any(r['valid'] for r in rows if r['run']==run and r['condition']==condition):
        ax.text(.5,.5,'No accepted output; NOT measured zero',ha='center',transform=ax.transAxes)
fig.suptitle('Offline comparison | separate same-setting video reference, NOT synchronized GT')
fig.tight_layout();fig.savefig(OUT/'comparison.png',dpi=140);plt.close(fig)

lines=['# 电钻球心与窗口离线对照结果','',
       '完成14个固定结束时刻 × 3种窗口 × 3种球心方案，共126条状态。生产算法未修改。', '',
       '## 本轮结论', '',
       '20ms + median3保留为下一轮候选。相对同档位名义参考，8个旋转片段中7个可用，平均绝对参考偏差11.52 RPM；10ms + fixed基线同7个可用片段为18.47 RPM。约38%的下降仅针对这组近似参考与已观察片段，不是独立真值精度提升认证。',
       'evaluation块候选与基线各4/4可用，平均绝对参考偏差12.21与21.92 RPM；development块各3/4，为10.60与13.85 RPM。此处是结果分层，不声称提前锁定候选后盲测。',
       '固定中心的20/40ms在3个有几何的纯平移窗口中均产生非零输出（37–67 RPM）；linear/median3在相同窗口没有接受非零输出。另1个平移窗几何失败，不能计作成功拒绝。静止两窗全部事件不足，不计作测准零。',
       '40ms + median3仅6/8旋转窗口可用，参考偏差12.01 RPM，未体现整体优势；不直接把窗口越长当作更好。本轮仅优化既有几何观测的时间处理，未修复缺测的轮廓前端。', '',
       '## 近似参考', '',
       '用户已接受60fps手机标注作为近似转速参考。独立从marker CSV重算确认视频1为249.4359 RPM、视频2为244.3868 RPM。',
       '用户本轮确认：另拍的同档位参考，转速近似相同。两段RPM平均246.9114作为名义档位参考，读数范围并集约233.88–261.53 RPM。该范围不是电钻跨录制漂移的标定界限。只报告相对档位参考偏差，不称同步真值误差；未把参考用于搜索、几何选择或调参。', '',
       '## 分组结果', '', '|录制/条件|窗口ms|球心|有效/预定|RPM中位数|有效范围|平均绝对参考偏差RPM|', '|---|---:|---|---:|---:|---|---:|']
for g in groups:
    med='—' if g['median'] is None else f"{g['median']:.1f}"
    extent='—' if g['minimum'] is None else f"{g['minimum']:.1f}–{g['maximum']:.1f}"
    deviation=g['approximate_reference_abs_deviation_mean_rpm']
    dev='—' if deviation is None else f'{deviation:.1f}'
    lines.append(f"|run{g['run']} {g['condition']}|{g['window_ms']}|{g['method']}|{g['valid']}/{g['planned']}|{med}|{extent}|{dev}|")
timings=[r['elapsed_s'] for r in rows if r['elapsed_s'] is not None]
lines += ['', '## 计算与解释边界', '',
          f"单次有输入估计耗时中位数{np.median(timings):.3f}s，最大{max(timings):.3f}s；包含球面提升和搜索，不含文件读取/原轮廓提取。",
          'fixed：窗口中点固定中心；linear：逐事件插值；median3：连续三点中值平滑后插值。三者同窗共享事件ID，40ms等长窗会因几何缺口失败，失败保留。',
          '各长度采用相同窗口结束时刻，但事件数和5000上限造成采样密度差异；结果是当前计算预算下整体表现，不能单独归因于时长。',
          '几何采用离线观测，可能含未来支点；lookahead字段是三方案共同的保守上界，非在线延迟。半径保持60px，未优化半径、外轮廓检测器或CMAX门限。',
          'development/evaluation按既有时间块区分，本次方案预先固定；历史已观察过这些片段，不能冒称独立盲测。',
          '未输出不等于测到零；转速稳定不等于轴正确。三维轴仅保留候选，未报告精确轴误差。', '',
          '## 复现', '',
          '`g1_spin/.venv/bin/python g1_spin/experiments/drill_window_geometry_20261110/run.py`',
          '首次运行输出拒绝覆盖。后续执行同目录 verify.py 与 summarize.py 核验并生成汇总。',
          'run.py保存运行前未配对的参考状态；reference.json另记录运行期间用户确认的同档位关系及参考汇总。重新开展实验时需显式保留该确认，summarize.py要求确认后的名义参考字段。',
          '来源与参数见protocol.md/freeze.json，逐条结果见results.csv/results.jsonl，独立输入核验见independent_verification.json。']
(OUT/'report.md').write_text('\n'.join(lines)+'\n')
print('Report, figure, per-run summaries, independent video reference check saved')
