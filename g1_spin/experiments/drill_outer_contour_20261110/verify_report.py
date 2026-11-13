"""Verify actual pixel support and raw spin inputs; report candidate limitations."""
import importlib.util
import json
from pathlib import Path
import numpy as np
OUT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('front',OUT/'run.py');f=importlib.util.module_from_spec(s);s.loader.exec_module(f)
rows=json.loads((OUT/'measurements.json').read_text());anchors=json.loads((f.SOURCE/'anchors.json').read_text());assert len(rows)==144
checks=[]
for r in rows:
    path=f.SOURCE/r['input_file'];assert f.old.digest(path)==r['input_sha256'];raw=np.load(path)['events']
    observed=set(zip(raw['x'].tolist(),raw['y'].tolist()));points=np.array(r['candidate']['points']).reshape(-1,2)
    assert all(tuple(p) in observed for p in points)
    assert len(set(r['candidate']['bins']))==len(points)
    a=next(a for a in anchors if a['case']==r['case']);anchor=np.array(a['center']);direction=np.array(a['direction']);perp=np.array([-direction[1],direction[0]])
    side=(points-anchor)@perp
    for key,mask in [('all',np.ones(len(points),bool)),('side_a',side<0),('side_b',side>=0)]:
        fit=r['fits'][key]
        if 'center' in fit:
            residual=np.sqrt(np.sum((points[mask]-fit['center'])**2,axis=1))-a['radius']
            assert np.isclose(np.sqrt(np.mean(residual**2)),fit['rms_px'])
    checks.append(dict(case=r['case'],mid_us=r['mid_us'],duration_us=r['duration_us'],actual_pixel_support_verified=True))
spin=json.loads((OUT/'spin_results.json').read_text());assert len(spin)==12
inputchecks=[]
for run in [3,5]:
    p2=f.base.p2;old=f.old;folder=p2.ARCHIVE/'decoded'/p2.NAMES[run];m=json.loads((folder/'manifest.json').read_text());parts=[np.load(folder/x,mmap_mode='r') for x in m['parts']];offsets=np.cumsum([0]+[len(p) for p in parts]);origin=int(parts[0]['t'][0])
    ots,oxy=p2.observations(run)
    def smooth(ts,xy):
        xy=xy.copy();result=xy.copy()
        for i in range(1,len(ts)-1):
            if np.all(np.diff(ts[i-1:i+2])<=30000):result[i]=np.sort(xy[i-1:i+2],axis=0)[1]
        return result
    oxy=smooth(ots,oxy)
    for a in [a for a in anchors if a['run']==run]:
        path=OUT/'spin_inputs'/f"{a['case']}.npz"
        if not path.exists():continue
        data=np.load(path);begin=a['end_us']-20000;end=a['end_us'];raw=[];ids=[]
        for i,p in enumerate(parts):
            if p['t'][-1]<origin+begin or p['t'][0]>=origin+end:continue
            mask=(p['t']>=origin+begin)&(p['t']<origin+end);raw.append(p[mask]);ids.append(np.flatnonzero(mask)+offsets[i])
        raw=np.concatenate(raw);ids=np.concatenate(ids);t=raw['t'].astype(np.int64)-origin;pix=np.column_stack([raw['x'],raw['y']]).astype(float)
        rr=[x for x in rows if x['case']==a['case'] and x['duration_us']==5000 and x['fits']['all']['valid']]
        nts=np.array([x['mid_us'] for x in rr]);nxy=smooth(nts,np.array([x['fits']['all']['center'] for x in rr]))
        cs={key:np.column_stack([np.interp(t,ts,xy[:,dim]) for dim in range(2)]) for key,ts,xy in [('old_median3',ots,oxy),('outer_median3',nts,nxy)]}
        mask=np.logical_and.reduce([np.sum((pix-c)**2,axis=1)<=(.8*a['radius'])**2 for c in cs.values()]);ix=np.flatnonzero(mask)
        assert np.array_equal(ids[ix],data['common_ids'])
        if len(ix)>5000:ix=ix[np.linspace(0,len(ix)-1,5000).astype(int)]
        assert np.array_equal(ids[ix],data['ids']) and np.array_equal(raw[ix],data['events'])
        for key in cs:assert np.allclose(cs[key][ix],data[key],rtol=0,atol=1e-9)
        inputchecks.append(a['case'])
    del parts
for fn in ['freeze.json','spin_freeze.json']:
    for path,h in json.loads((OUT/fn).read_text()).items():assert f.old.digest(f.base.ROOT/path)==h

def good(x,key):return x[key]['side_a']['valid'] and x[key]['side_b']['valid']
def gap(x,key):return float(np.linalg.norm(np.array(x[key]['side_a']['center'])-x[key]['side_b']['center']))
def med(v):return float(np.median(v)) if v else None
summary=[]
for a in anchors:
    for duration in [1000,5000]:
        rr=[x for x in rows if x['case']==a['case'] and x['duration_us']==duration];paired=[x for x in rr if good(x,'fits') and good(x,'old_fits')]
        for key in ['old_fits','fits']:
            summary.append(dict(case=a['case'],duration_us=duration,method=key,paired_sides=sum(good(x,key) for x in rr),
                 all_valid=sum(x[key]['all']['valid'] for x in rr),common_count=len(paired),common_side_gap_median_px=med([gap(x,key) for x in paired])))
prior=json.loads((OUT.parent/'drill_geometry_diagnosis_20261110/geometry_audit.json').read_text());spin_summary=[]
for r in spin:
    res=r['result'];angle=None
    if r['run']==3 and res['valid']:
        p=next(x for x in prior if x['run']==3 and x['end_us']==int(r['case'].split('_')[1]))['shaft']['full20']
        if p['valid']:
            aa=np.array(res['axis_camera'][:2]);bb=np.array(p['direction']);angle=float(np.degrees(np.arccos(np.clip(abs(aa@bb)/np.linalg.norm(aa)/np.linalg.norm(bb),-1,1))))
    spin_summary.append(dict(case=r['case'],run=r['run'],method=r['method'],valid=res['valid'],rpm=res['spin_rpm'],shaft_projection_deg=angle,reason=res['failure_reason']))
f.old.save(OUT/'summary.json',dict(geometry=summary,spin=spin_summary))
f.old.save(OUT/'verification.json',dict(exposures=144,actual_pixels_verified=144,spin_rows=12,raw_spin_inputs_verified=inputchecks,source_hashes_unchanged=True))
print(json.dumps(dict(spin=spin_summary,geometry_summary=summary),indent=2))
# Supplement independently checks raw physical IDs and window membership.
standalone=json.loads((OUT/'old_standalone.json').read_text());data=np.load(OUT/'spin_inputs/old_standalone.npz')
p2=f.base.p2;folder=p2.ARCHIVE/'decoded'/p2.NAMES[3];manifest=json.loads((folder/'manifest.json').read_text());offset=0;recovered=[]
for name in manifest['parts']:
    part=np.load(folder/name,mmap_mode='r');sel=data['ids'][(data['ids']>=offset)&(data['ids']<offset+len(part))]-offset
    recovered.append(part[sel]);offset+=len(part)
assert np.array_equal(np.concatenate(recovered),data['events'])
times=data['events']['t'].astype(np.int64)-int(data['origin_us']);assert np.all((times>=4315000)&(times<4335000))
ts,xy=p2.observations(3);xy=smooth(ts,xy);centers=np.column_stack([np.interp(times,ts,xy[:,i]) for i in range(2)])
assert np.allclose(centers,data['centers'],rtol=0,atol=1e-9)
verified=json.loads((OUT/'verification.json').read_text());verified['standalone_old_raw_ids_and_centers_verified']=True
f.old.save(OUT/'verification.json',verified)

slant=[]
for dur in [1000,5000]:
    rr=[x for x in rows if x['run']==3 and x['duration_us']==dur];paired=[x for x in rr if good(x,'fits') and good(x,'old_fits')]
    for k in ['old_fits','fits']:
        slant.append(dict(duration_us=dur,method=k,coverage=sum(good(x,k) for x in rr),planned=len(rr),common=len(paired),gap=med([gap(x,k) for x in paired])))
lines=['# 外侧支撑弧：几何改善，尚不替换默认','',
       '完成144个快门的固定候选对照、6个20ms窗口的12条旋转/平移状态，以及1条独立旧前端补测。全程没有按参考转速或轴结果调整参数。', '',
       '## 实现', '',
       '每个角区选择最外侧、具有局部唯一像素和相邻角区支持的径向带；剔除孤立热点与固定ROI边界候选，保留实际观测像素，不向理想圆投影。随后沿用原固定半径拟合和质量检查。较旧方案同时改变了r±8到r±12搜索带，结果属于整个候选前端效果，不单独归因于外侧选择。', '',
       '## 几何结果', '',
       '|斜轴快门|前端|两侧同时有效/36|共同有效时刻数|共同球心差中位px|', '|---|---|---:|---:|---:|']
for x in slant:lines.append(f"|{x['duration_us']/1000:g}ms|{x['method']}|{x['coverage']}/36|{x['common']}|{x['gap']:.2f}|")
lines += ['', '5ms有明显改善；1ms支撑稀疏时该方法不能同等适用。两侧更一致仍可能存在共同偏差，未认证球心真值。', '',
          '## 旋转验证（同窗同事件ID）', '',
          '|片段|前端|状态|RPM|二维杆投影差°|','|---|---|---|---:|---:|']
for x in spin_summary:
    rpm='—' if x['rpm'] is None else f"{x['rpm']:.2f}";ang='—' if x['shaft_projection_deg'] is None else f"{x['shaft_projection_deg']:.2f}"
    lines.append(f"|{x['case']}|{x['method']}|{'有效' if x['valid'] else x['reason']}|{rpm}|{ang}|")
lines += ['',
          '三个共同有效斜轴窗的二维杆投影差中位：旧1.92°，新0.51°；转速中位：旧261.96RPM，新265.29RPM，没有同时改善。run3该方向没有单独确认的同步RPM真值，不用约247RPM选胜者。',
          f"run3_04335000新轮廓在窗口起点之前没有有效全环中心，导致共同准备失败；独立补测旧前端仍有效，{standalone['result']['spin_rpm']:.2f}RPM。因此正确覆盖是旧4/4、新3/4，不能把配对记录双方失败解释为可用率相同。",
          '两个纯平移窗两法均返回no_rotation_evidence，无新增非零输出。几何两侧差在第一个平移上下文有增大，说明分区一致性改善不普遍。样本仅两个，不外推完整误报率。', '',
          '## 决定与边界', '',
          '保留outer_supported为离线实验候选，不替换生产。它对5ms斜轴纹理干扰显示改进，但在稀疏1ms或部分边缘弱的快门上丢失支撑；外侧候选也可能选到外部反光或背景，不能保证silhouette身份。',
          '缺失几何不填、不自动切到结果更好的侧、不按已知轴纠正。保留旧前端和此前几何缺口回退；若后续整合，应先明确新前端几何不足时的回退与诊断，不将低残差等同可靠中心。',
          '下一步最有价值的是在新的录制/未用于本轮的上下文检验同一固定前端，并分析丢失的4.335s之前轮廓支持；当前不应反复在四个斜轴窗上调参数。', '',
          '## 验证和复现', '',
          '3个解析/热点单元测试通过；144快门选点均为原始观测像素、拟合残差重算通过；5个共同可计算窗口的原始ID、公共掩膜与新旧中心独立重建通过，独立旧前端补测输入另行核验；冻结来源哈希未变。',
          '使用g1_spin/.venv/bin/python：先pytest test_contour.py，再run.py、spin.py、supplement.py、verify_report.py。实验输出拒绝覆盖。',
          'contour.py为候选实现；measurements.json保存几何；spin_results.json与old_standalone.json保存旋转/失败状态；verification.json为核验；*_review.png为原始事件与实际选点。']
(OUT/'report.md').write_text('\n'.join(lines)+'\n')
print('Final report and standalone verification saved')
