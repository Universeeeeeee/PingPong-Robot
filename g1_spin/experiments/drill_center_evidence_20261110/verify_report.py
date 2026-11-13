"""Raw/geometry verification and fixed-evidence summary, no spin tuning."""
import importlib.util
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('evidence_run',OUT/'run.py');r=importlib.util.module_from_spec(s);s.loader.exec_module(r)
rows=json.loads((OUT/'measurements.json').read_text());anchors=json.loads((OUT/'anchors.json').read_text());assert len(rows)==144
assert len({(x['case'],x['mid_us'],x['duration_us']) for x in rows})==144
checks=[]
for run in r.CASES:
    folder=r.p2.ARCHIVE/'decoded'/r.p2.NAMES[run];manifest=json.loads((folder/'manifest.json').read_text())
    parts=[np.load(folder/f,mmap_mode='r') for f in manifest['parts']];offsets=np.cumsum([0]+[len(p) for p in parts]);origin=int(parts[0]['t'][0])
    for row in [x for x in rows if x['run']==run]:
        a=next(x for x in anchors if x['case']==row['case']);center=np.array(a['center']);radius=a['radius']
        path=OUT/row['input_file'];assert r.old.digest(path)==row['input_sha256'];data=np.load(path)
        allids=[];allraw=[];begin=row['mid_us']-row['duration_us']//2;end=row['mid_us']+row['duration_us']//2
        for i,p in enumerate(parts):
            if p['t'][-1]<origin+begin or p['t'][0]>=origin+end:continue
            mask=(p['t']>=origin+begin)&(p['t']<origin+end)
            q=p[mask];ids=np.flatnonzero(mask)+offsets[i]
            delta=np.column_stack([q['x'],q['y']]).astype(float)-center
            local=(delta**2).sum(axis=1)<=(radius+190)**2
            allraw.append(q[local]);allids.append(ids[local])
        assert np.array_equal(np.concatenate(allids),data['ids']) and np.array_equal(np.concatenate(allraw),data['events'])
        pts=np.array(row['ring_points']).reshape(-1,2);d=np.array(a['direction']);perp=np.array([-d[1],d[0]])
        sides=(pts-center)@perp
        for method,mask in [('all',np.ones(len(pts),bool)),('side_a',sides<0),('side_b',sides>=0)]:
            fit=row['fits'][method]
            if 'center' in fit:
                res=np.sqrt(np.sum((pts[mask]-fit['center'])**2,axis=1))-radius
                assert np.isclose(np.sqrt(np.mean(res**2)),fit['rms_px'])
        timing=row['bin_median_time_offset_us'];assert all(abs(t)<=row['duration_us']/2 for t in timing)
        checks.append(dict(case=row['case'],mid_us=row['mid_us'],duration_us=row['duration_us'],raw_ids_verified=True))
    del parts
# Analytic independent circle fixture exercises fit, not the DVS observation model.
theta=np.linspace(0,2*np.pi,72,endpoint=False);true=np.array([103.,98.]);points=true+60*np.column_stack([np.cos(theta),np.sin(theta)])
errors=[]
for subset in [points,points[points[:,0]<true[0]],points[points[:,0]>=true[0]]]:
    f=r.fit_circle(subset,np.array([100.,100.]),60,6);assert f['valid'];errors.append(float(np.linalg.norm(np.array(f['center'])-true)))
assert max(errors)<1e-5
for filename,h in json.loads((OUT/'freeze.json').read_text()).items():assert r.old.digest(r.ROOT/filename)==h
r.old.save(OUT/'verification.json',dict(exposures=144,raw_verified=len(checks),circle_fixture_max_error_px=max(errors),
                                    frozen_sources_unchanged=True,checks=checks))
summaries=[];fig,axes=plt.subplots(4,2,figsize=(13,12))
for a,ax in zip(anchors,axes.flat):
    center=np.array(a['center']);d=np.array(a['direction']);perp=np.array([-d[1],d[0]])
    for duration in [1000,5000]:
        rr=[x for x in rows if x['case']==a['case'] and x['duration_us']==duration]
        pair=[x for x in rr if all(x['fits'][k]['valid'] for k in ['side_a','side_b'])]
        gaps=[float(np.linalg.norm(np.array(x['fits']['side_a']['center'])-x['fits']['side_b']['center'])) for x in pair]
        common=[x for x in pair if x['shaft']['valid'] and x['old_center'] is not None]
        traces={}
        for key in ['side_a','side_b','old','shaft']:
            values=[]
            for x in common:
                c=x['fits'][key]['center'] if key in ['side_a','side_b'] else x['old_center'] if key=='old' else None
                values.append(float(np.dot(np.array(c)-center,perp)) if c is not None else x['shaft']['offset_px'])
            traces[key]=np.array(values)
        metrics={}
        if len(common)>=3:
            target=traces['shaft']-traces['shaft'].mean()
            for key in ['side_a','side_b','old']:
                v=traces[key];centered=v-v.mean()
                metrics[key]=dict(span_px=float(np.ptp(v)),shaft_span_px=float(np.ptp(traces['shaft'])),
                    demeaned_difference_rms_px=float(np.sqrt(np.mean((centered-target)**2))),
                    correlation=float(np.corrcoef(v,traces['shaft'])[0,1]) if np.std(v)>1e-9 and np.std(traces['shaft'])>1e-9 else None)
        summaries.append(dict(case=a['case'],duration_us=duration,planned=len(rr),
            side_a_valid=sum(x['fits']['side_a']['valid'] for x in rr),side_b_valid=sum(x['fits']['side_b']['valid'] for x in rr),
            shaft_valid=sum(x['shaft']['valid'] for x in rr),paired_sides=len(pair),
            side_center_distance_median_px=float(np.median(gaps)) if gaps else None,
            common_count=len(common),common_times_us=[x['mid_us'] for x in common],transverse_metrics=metrics))
        if duration==5000:
            # Plot each signal on its own available support; don't fill missing values.
            for key,color in [('old','#333333'),('side_a','#20a040'),('side_b','#bd40ad'),('shaft','#2089bc')]:
                times=[];values=[]
                for x in rr:
                    if key=='shaft':value=x['shaft'].get('offset_px') if x['shaft']['valid'] else None
                    else:
                        c=x['old_center'] if key=='old' else x['fits'][key].get('center') if x['fits'][key]['valid'] else None
                        value=float((np.array(c)-center)@perp) if c is not None else None
                    times.append((x['mid_us']-a['end_us'])/1000);values.append(np.nan if value is None else value)
                v=np.array(values);v-=np.nanmean(v) if np.isfinite(v).any() else 0
                ax.plot(times,v,'o-',ms=3,color=color,label=key)
            ax.set(title=a['case'],xlabel='Relative to window end (ms)',ylabel='Demeaned transverse offset (px)');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('5ms raw geometry: separate support; demeaned traces are not ground truth')
fig.tight_layout();fig.savefig(OUT/'transverse_comparison.png',dpi=130);plt.close(fig)
r.old.save(OUT/'summary.json',summaries)
lines=['# 球心变化的原始事件交叉检查','',
       '完成8个上下文、每个9时刻×1/5ms，共144个原始事件快门。没有运行旋转优化，也没有修改生产算法。', '',
       '## 判断', '',
       '证据支持旧球心轨迹中包含轮廓拟合干扰；不能把全部变化解释成真实振动，也不能宣布全部是假运动。尤其斜轴段，两侧轮廓与外杆没有形成完全一致的中心轨迹。当前不宜直接把这条中心轨迹求速度并用于强平移补偿。',
       '四个斜轴上下文在5ms共同支持时刻中，旧球心横向范围6.25–7.77px，杆横向范围1.17–2.33px；旧中心相对杆的去均值差RMS为1.75–2.40px。比较只约束杆垂直方向，不能测沿杆运动、杆倾斜或球安装偏心。',
       '两侧轮廓的中心距离中位约3.08–4.45px（斜轴5ms），尽管单侧拟合RMS通过。低圆拟合残差不等于真实球心正确。原始图可见旋转纹理靠近候选外环，存在混入风险，但未给每条边缘做外轮廓身份标注。',
       '纯平移run5_04135000中，两侧与杆横向变化的相关系数约0.99–1.00，去均值差RMS约0.60–0.83px，说明这套辅助检查能在部分正例中看到一致运动；另一平移上下文两侧中心偏差更大，不据此泛化。',
       '1ms图更稀疏，5ms图有更多支撑但也更易混入随旋转变化的纹理；不把短快门自动视为真值。当前新候选由固定先验圈选，仍有选择偏差，不能以某一侧更贴近杆为理由自动选择该侧。', '',
       '## 5ms结果（全部失败留在分母）', '',
       '|上下文|两侧同时有效/9|两侧球心距离中位px|与杆及旧中心共同支持数|旧中心横向范围px|杆横向范围px|',
       '|---|---:|---:|---:|---:|---:|']
for x in summaries:
    if x['duration_us']!=5000:continue
    v=x['transverse_metrics'].get('old');fmt=lambda z:'—' if z is None else f'{z:.2f}'
    lines.append(f"|{x['case']}|{x['paired_sides']}/9|{fmt(x['side_center_distance_median_px'])}|{x['common_count']}|{fmt(v['span_px'] if v else None)}|{fmt(v['shaft_span_px'] if v else None)}|")
lines += ['', '范围和RMS只在同一行的共同有效时刻计算，不跨不同支持集合比较。相关系数来自4–8个时刻，属于描述性诊断，不是统计精度认证。', '',
          '## 下一步实施方向', '',
          '先构建与旋转纹理解耦的轮廓前端：明确外轮廓支持，并以分区一致性和可辨几何约束评估中心。现有两侧检查用于暴露问题，不直接设成自动质量门槛或替换中心。若需要定量判定真实振动幅度，须同次同步拍摄外部球心参考；现有另拍60fps只提供档位转速。',
          '保留已验证的几何缺口回退，保留固定中心与补偿候选两路；本轮没有足够证据推广新的默认球心算法。', '',
          '## 核验与产物', '',
          f'144快门的原始事件ID、时间范围、固定ROI与原官方分块重建一致；拟合残差独立复算通过。解析整圆/半圆已知中心夹具最大误差{max(errors):.3g}px，仅验证求解数学。冻结文件哈希未变。',
          'measurements.json为逐快门原始测量，summary.json为支持计数和横向一致性；transverse_comparison.png为位移曲线；各*_raw.png为原始1/5ms事件与候选圆。raw/保留官方事件子集及物理ID。',
          '复现：使用g1_spin/.venv/bin/python执行本目录run.py（拒绝覆盖），再执行verify_report.py。']
(OUT/'report.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(dict(exposures=len(rows),verified=len(checks),circle_fixture_max_error_px=max(errors),report=str(OUT/'report.md'))))
