"""Rebuild raw masks, check geometry-only route, report bounded diagnosis."""
import importlib.util
import json
from pathlib import Path
from collections import Counter
import numpy as np
from route import covered

OUT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('diagnosis',OUT/'run.py');d=importlib.util.module_from_spec(s);s.loader.exec_module(d)
rows=[json.loads(line) for line in (OUT/'results.jsonl').read_text().splitlines()]
assert len(rows)==48 and len({(r['window'],r['method']) for r in rows})==48
checks=[]
for run in d.ENDS:
    folder=d.p2.ARCHIVE/'decoded'/d.p2.NAMES[run];m=json.loads((folder/'manifest.json').read_text())
    parts=[np.load(folder/f,mmap_mode='r') for f in m['parts']];offsets=np.cumsum([0]+[len(x) for x in parts]);origin=int(parts[0]['t'][0])
    ts,xy=d.p2.observations(run);sm=xy.copy()
    for i in range(1,len(ts)-1):
        if max(np.diff(ts[i-1:i+2]))<=30000:sm[i]=np.sort(xy[i-1:i+2],axis=0)[1]
    for path in sorted((OUT/'inputs').glob(f'run{run}_*.npz')):
        data=np.load(path);group=[r for r in rows if r['window']==path.stem];assert len(group)==3
        r=group[0];end=r['end_us'];begin=end-r['window_ms']*1000
        raw=[];rawids=[]
        for i,part in enumerate(parts):
            mask=(part['t']>=origin+begin)&(part['t']<origin+end)
            raw.append(part[mask]);rawids.append(np.flatnonzero(mask)+offsets[i])
        raw=np.concatenate(raw);rawids=np.concatenate(rawids);times=raw['t'].astype(np.int64)-origin
        pixels=np.column_stack([raw['x'],raw['y']]);radius=float(data['radius'])
        fixed=np.array([np.interp((begin+end)//2,ts,xy[:,dim]) for dim in range(2)])
        centers={'fixed':np.tile(fixed,(len(raw),1)),
                 'linear':np.column_stack([np.interp(times,ts,xy[:,dim]) for dim in range(2)]),
                 'median3':np.column_stack([np.interp(times,ts,sm[:,dim]) for dim in range(2)])}
        keep=np.ones(len(raw),bool)
        for center in centers.values():keep&=np.sum((pixels-center)**2,axis=1)<=(.8*radius)**2
        ix=np.flatnonzero(keep);assert np.array_equal(rawids[ix],data['common_ids'])
        if len(ix)>5000:ix=ix[np.linspace(0,len(ix)-1,5000).astype(int)]
        assert np.array_equal(rawids[ix],data['ids']) and np.array_equal(raw[ix],data['events'])
        for method in centers:assert np.allclose(centers[method][ix],data[method],rtol=0,atol=1e-9)
        for item in group:
            result=item['result']
            if result['valid']:assert np.isclose(np.linalg.norm(result['omega_rad_s'])*60/(2*np.pi),result['spin_rpm'])
            else:assert result['spin_rpm'] is None
        checks.append(dict(window=path.stem,events=len(ix),raw_common_mask_geometry_verified=True))
    del parts

routes=json.loads((OUT/'routing.json').read_text());assert len(routes)==44
for r in routes:
    if r['condition']=='static':continue
    ts,_=d.p2.observations(r['run']);end=r['end_us']
    def independent_coverage(begin):
        if begin<ts[0] or end>ts[-1]:return False
        return all(not (b-a>30000 and a<end and b>begin) for a,b in zip(ts,ts[1:]))
    expected='20ms_median3' if independent_coverage(end-20000) else '10ms_fixed_fallback' if independent_coverage(end-10000) else 'geometry_unavailable'
    assert r['route']==expected
    if expected!='geometry_unavailable':assert r['reused_result'] is not None
for filename,h in json.loads((OUT/'freeze.json').read_text()).items():assert d.old.digest(d.ROOT/filename)==h
verification=dict(result_rows=48,prepared_windows=len(checks),routes=44,route_tests=5,checks=checks,protected_hashes_unchanged=True)
d.old.save(OUT/'verification.json',verification)

def angle(a,b):
    return float(np.degrees(np.arctan2(abs(np.cross(np.asarray(a),np.asarray(b))),abs(np.dot(a,b)))))
audit=json.loads((OUT/'geometry_audit.json').read_text());table=[]
for ms in [10,20]:
    for method in d.METHODS:
        group=[r for r in rows if r['run']==3 and r['window_ms']==ms and r['method']==method]
        good=[r for r in group if r['result']['valid']];angles=[]
        for r in good:
            shaft=next(a['shaft']['full20'] for a in audit if a['run']==3 and a['end_us']==r['end_us'])
            if shaft['valid']:angles.append(angle(r['result']['axis_camera'][:2],shaft['direction']))
        table.append(dict(window_ms=ms,method=method,valid=len(good),planned=len(group),
                          rpm_median=float(np.median([r['result']['spin_rpm'] for r in good])),
                          shaft_projection_median_deg=float(np.median(angles)),angles_deg=angles))
shaft_checks=[]
for a in audit:
    if a['run']==3:
        first,last=a['shaft']['first10'],a['shaft']['last10']
        shaft_checks.append(dict(end_us=a['end_us'],first_last_angle_deg=angle(first['direction'],last['direction']) if first['valid'] and last['valid'] else None))
route_summary=[]
nominal=json.loads((d.prev.old.OUT/'reference.json').read_text())['nominal_same_setting_rpm']
for run in [1,2,3,4,5]:
    for condition in (['translation','static'] if run==5 else ['rotation']):
        rr=[r for r in routes if r['run']==run and r['condition']==condition]
        good=[r for r in rr if r['reused_result'] and r['reused_result']['valid']]
        route_summary.append(dict(run=run,condition=condition,planned=len(rr),valid=len(good),
             routes=dict(Counter(r['route'] for r in rr)),reference_abs_deviation_mean_rpm=float(np.mean([abs(r['reused_result']['spin_rpm']-nominal) for r in good])) if good and run in [1,4] else None))
d.old.save(OUT/'summary.json',dict(factorial=table,shaft_reference_temporal_checks=shaft_checks,routing=route_summary,
    decision='Geometry fallback is reproducible; tilted-axis issue remains; do not promote default'))
lines=['# 几何缺口与斜轴诊断结果','',
       '完成8个既有问题/控制区间 × 2种窗口 × 3种中心，共48条状态；另在前轮44区间上验证只依赖几何的路由。固定范围诊断，不是新独立采集。', '',
       '## 结论', '',
       '1. 两个漏窗来自40ms几何观测空洞。20ms不覆盖时回退到有完整几何的10ms固定中心，可以恢复两窗；不填缺口、不放宽30ms限制、不看旋转得分决定回退。',
       '2. 斜轴的补偿几何是优先排查对象。10ms下fixed/linear/median3杆投影差中位约0.24/5.08/2.31度；20ms约1.18/2.27/2.65度。固定时长也能出现差异，不能只归因延长窗口；但无法据此唯一判定真实中心轨迹谁更准。',
       '3. 前后10ms杆投影差约0.25–1.29度；弱参考本身有波动，但这四窗没有显示足以解释所有中心方案差异的统一漂移。杆投影来自同传感器且不是三维轴真值。',
       '4. 不能统一退回固定中心：两个纯平移控制在20ms fixed误报46.85/66.41RPM，而linear/median3未接受非零。需要兼顾平移补偿与轴偏差。',
       '5. 保留实验候选及几何回退，不推广生产默认。后续应优先独立检查中心轨迹/有效观测时刻，不能通过约247RPM择峰或按钻杆方向强行修轴。', '',
       '## 斜轴六组对照', '', '|窗口ms|中心|有效/预定|RPM中位|杆投影差中位°|', '|---:|---|---:|---:|---:|']
for r in table:lines.append(f"|{r['window_ms']}|{r['method']}|{r['valid']}/{r['planned']}|{r['rpm_median']:.2f}|{r['shaft_projection_median_deg']:.2f}|")
lines += ['', '同一时长三种中心共享事件ID；跨时长不同样本且有5000上限，不是纯时长效应。与上一轮差异还含共同掩膜变化，未覆盖旧结果。run3没有独立三维轴和该录制的同步RPM标签。', '',
          '## 几何路由回放', '', '|run/条件|有效/预定|路由计数|', '|---|---:|---|']
for r in route_summary:lines.append(f"|{r['run']}/{r['condition']}|{r['valid']}/{r['planned']}|{r['routes']}|")
lines += ['',
          '旋转24/32可输出（原候选22/32，原基线经补测22/32）；纯平移0/8接受非零；静止0/4因事件不足，不能视作准确测零。所有数字来自记录级复用，不是重新执行完整在线流水线。',
          '回退两个窗口使用前轮独立10ms基线的实际结果249.18/237.95RPM；完整20ms仍使用前轮候选，斜轴问题并未被路由修复。严格检查区间内部缺口，避免只检首尾产生自洽假象。',
          '几何支点含离线未来观测；实时性未验证。路由不会在CMAX不收敛/拒绝之后二次挑结果。', '',
          '## 核验与复现', '',
          '5项路由边界测试通过；48状态/14份输入的原始事件ID、公共掩膜、中心插值与平滑独立重建通过；44项路由由逐间隔检查交叉核对；冻结源码及旧实验哈希未变。',
          '`g1_spin/.venv/bin/python -m pytest -q g1_spin/experiments/drill_geometry_diagnosis_20261110/test_route.py`',
          '`g1_spin/.venv/bin/python g1_spin/experiments/drill_geometry_diagnosis_20261110/run.py`（拒绝覆盖）',
          '`g1_spin/.venv/bin/python g1_spin/experiments/drill_geometry_diagnosis_20261110/verify_report.py`',
          'source: protocol.md、freeze.json、geometry_audit.json、routing.json、results.jsonl、verification.json。']
(OUT/'report.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(dict(verified_inputs=len(checks),factorial=table,shaft=shaft_checks,routing=route_summary),indent=2))
