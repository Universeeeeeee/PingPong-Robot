"""Unmodified strict local-plane flow on complete synthetic ball events."""
from pathlib import Path
import importlib.util
import json
import sys
from time import perf_counter

import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[3]
PHASE3 = ROOT / 'g1_spin/experiments/drill_phase3_20261106'
sys.path.insert(0, str(ROOT / 'g1_spin/src'))
from g1_spin.events import EventBatch
from g1_spin.optical_flow import estimate_flow, solve_constraints

specification = importlib.util.spec_from_file_location('strict_stage3_flow', PHASE3 / 'flow.py')
flow = importlib.util.module_from_spec(specification)
specification.loader.exec_module(flow)
protocol = json.loads((OUT / 'protocol.json').read_text())
results = []
for scene in protocol['scenes']:
    name = f"rpm{scene['rpm']}_ms{int(scene['duration_ms'])}"
    data = np.load(OUT / f'{name}_events.npz')
    truth = np.load(OUT / f'{name}_truth.npz')['omega']
    xy = data['xy']
    times = data['times_us'] * 1e-6
    started = perf_counter()
    planes, meta = flow.extract(xy, times, data['polarity'], np.array(protocol['center']), protocol['radius'])
    elapsed = perf_counter()-started
    initial = flow.initialize(planes['points'], planes['gradients'], protocol['radius'], 'S')
    native = solve_constraints(planes['points'], planes['gradients'], protocol['radius'])
    batch = EventBatch(xy[:,0], xy[:,1], data['times_us'], data['polarity'])
    actual = estimate_flow(batch, *protocol['center'], protocol['radius'])
    if meta['accepted_planes'] >= 12:
        assert actual['valid'] == native['valid']
        assert actual['failure_reason'] == native['failure_reason']
        if actual['valid']:
            assert np.allclose(actual['omega_rad_s'], native['omega_rad_s'], atol=1e-10)
    else:
        assert actual['failure_reason'] == 'insufficient_stable_local_planes'
    row = dict(scene=name,events=len(xy),source='All generated sphere events, no training/holdout split or event cap; diagnostic only',
               metadata=meta,initial=initial,production=actual,extraction_s=elapsed)
    if initial['available']:
        omega = np.array(initial['raw_omega_rad_s'])
        row['candidate_rpm'] = float(np.linalg.norm(omega)*60/(2*np.pi))
        row['candidate_rpm_error_percent'] = float(100*(np.linalg.norm(omega)/np.linalg.norm(truth)-1))
        row['candidate_axis_error_deg'] = float(np.degrees(np.arccos(np.clip(omega@truth/(np.linalg.norm(omega)*np.linalg.norm(truth)),-1,1))))
    results.append(row)
    print(name,'events',len(xy),'planes',meta['accepted_planes'],'production',actual['valid'],actual['failure_reason'],
          'rpm',row.get('candidate_rpm'),'axis',row.get('candidate_axis_error_deg'),flush=True)
(OUT/'strict_flow_results.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
lines = ['## 原严格光流前端的补充诊断','',
         '对上述六场景完整球内事件调用原production estimate_flow，并与phase3 flow.extract＋S逐项核对。保持4px、±2ms、300seed、至少10邻居、空间比0.08、RMS≤max(20us,0.15std)、至少12条全局约束等全部原门限。此项使用全事件，属于前端能力诊断，未用来替代前表训练/留出的估计结果。','',
         '| 场景 | 总事件 | 接受平面/300 | production输出 | S候选RPM | S候选轴误差 |','|---|---:|---:|---|---:|---:|']
for row in results:
    rpm = f"{row['candidate_rpm']:.2f}" if 'candidate_rpm' in row else '—'
    angle = f"{row['candidate_axis_error_deg']:.3f}°" if 'candidate_axis_error_deg' in row else '—'
    status = '有效' if row['production']['valid'] else row['production']['failure_reason']
    lines.append(f"| {row['scene']} | {row['events']} | {row['metadata']['accepted_planes']} | {status} | {rpm} | {angle} |")
lines += ['', '全部场景production入口与stage3提取/求解结果一致。有效S候选由事件推导后再与生成真值比较，真值未进入前端或最小二乘。具体失败计数及残差见strict_flow_results.json。']
(OUT/'strict_flow_report.md').write_text('\n'.join(lines)+'\n')
