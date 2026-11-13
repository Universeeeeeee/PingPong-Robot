import importlib.util,json
from pathlib import Path
from collections import Counter
import numpy as np
OUT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('trial',OUT/'run.py');trial=importlib.util.module_from_spec(s);s.loader.exec_module(trial)
rows=json.loads((OUT/'results.json').read_text());geometries=json.loads((OUT/'geometry.json').read_text());routes=json.loads((OUT/'routes.json').read_text())
assert len(rows)==36 and len({(r['case'],r['lane']) for r in rows})==36
checks=[]
def smooth(ts,xy):
    result=xy.copy()
    for i in range(1,len(ts)-1):
        if np.all(np.diff(ts[i-1:i+2])<=30000):result[i]=np.sort(xy[i-1:i+2],axis=0)[1]
    return result
for run in trial.BLOCKS:
    folder=trial.p2.ARCHIVE/'decoded'/trial.p2.NAMES[run];m=json.loads((folder/'manifest.json').read_text());parts=[np.load(folder/x,mmap_mode='r') for x in m['parts']];offsets=np.cumsum([0]+[len(p) for p in parts]);origin=int(parts[0]['t'][0])
    def fetch(begin,end):
        chunks=[];ids=[]
        for i,p in enumerate(parts):
            if p['t'][-1]<origin+begin or p['t'][0]>=origin+end:continue
            mask=(p['t']>=origin+begin)&(p['t']<origin+end);chunks.append(p[mask]);ids.append(np.flatnonzero(mask)+offsets[i])
        return np.concatenate(chunks),np.concatenate(ids)
    for g in [g for g in geometries if g['run']==run]:
        path=OUT/g['input_file'];assert trial.old.digest(path)==g['input_sha256'];data=np.load(path)
        raw,ids=fetch(g['mid_us']-2500,g['mid_us']+2500);pix=np.column_stack([raw['x'],raw['y']]).astype(float)
        mask=np.sum((pix-g['anchor'])**2,axis=1)<=(g['radius']+190)**2
        assert np.array_equal(ids[mask],data['ids']) and np.array_equal(raw[mask],data['events'])
        observed=set(zip(data['events']['x'].tolist(),data['events']['y'].tolist()))
        assert all(tuple(p) in observed for p in g['candidate']['points'])
    ots,oxy=trial.p2.observations(run)
    for r in [r for r in rows if r['run']==run and r['input_file'] and not r['cached_same_route']]:
        path=OUT/r['input_file'];assert trial.old.digest(path)==r['input_sha256'];data=np.load(path);end=r['end_us'];begin=end-(10000 if r['route']=='old10' else 20000)
        raw,ids=fetch(begin,end);t=raw['t'].astype(np.int64)-origin;pix=np.column_stack([raw['x'],raw['y']]).astype(float)
        if r['route']=='old10':cs=np.tile([np.interp(end-5000,ots,oxy[:,d]) for d in range(2)],(len(t),1))
        elif r['route']=='old20' and r['condition']=='static':cs=np.tile([714.,355.],(len(t),1))
        else:
            ts,xy=ots,oxy
            if r['route']=='outer20':
                gg=[g for g in geometries if g['case']==r['case'] and g['fit']['valid']];ts=np.array([g['mid_us'] for g in gg]);xy=np.array([g['fit']['center'] for g in gg])
            xy=smooth(ts,xy);cs=np.column_stack([np.interp(t,ts,xy[:,d]) for d in range(2)])
        ix=np.flatnonzero(np.sum((pix-cs)**2,axis=1)<=(.8*trial.RADIUS[run])**2)
        assert np.array_equal(ids[ix],data['common_ids'])
        if len(ix)>5000:ix=ix[np.linspace(0,len(ix)-1,5000).astype(int)]
        assert np.array_equal(ids[ix],data['ids']) and np.array_equal(raw[ix],data['events']) and np.allclose(cs[ix],data['centers'],rtol=0,atol=1e-9)
        checks.append(r['input_file'])
    for r in [r for r in routes if r['run']==run]:
        nts=np.array([g['mid_us'] for g in geometries if g['case']==r['case'] and g['fit']['valid']]);end=r['end_us']
        def covered(ts,begin):
            return len(ts)>=2 and begin>=ts[0] and end<=ts[-1] and all(not (b-a>30000 and a<end and b>begin) for a,b in zip(ts,ts[1:]))
        old='old20' if r['condition']=='static' or covered(ots,end-20000) else 'old10' if covered(ots,end-10000) else 'unavailable'
        new='outer20' if covered(nts,end-20000) else old
        assert (old,new)==(r['old_route'],r['new_route'])
for r in rows:
    if r['cached_same_route']:
        other=next(x for x in rows if x['case']==r['case'] and x['lane']=='old');assert other['route']==r['route'] and other['result']==r['result']
for name,root in [('freeze.json',trial.base.ROOT),('routing_freeze.json',OUT)]:
    for p,h in json.loads((OUT/name).read_text()).items():assert trial.old.digest(root/p)==h
trial.old.save(OUT/'verification.json',dict(geometry_exposures=len(geometries),raw_spin_inputs=len(checks),states=36,routes_verified=18,tests=6,protected_sources_unchanged=True))

def angle(a,b,unsigned=False):
    a,b=np.asarray(a),np.asarray(b);c=float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)));return float(np.degrees(np.arccos(np.clip(abs(c) if unsigned else c,-1,1))))
def med(x):return float(np.median(x)) if x else None
nominal=json.loads((trial.old.OUT/'reference.json').read_text())['nominal_same_setting_rpm'];summary=[]
for run in trial.BLOCKS:
    for condition in (['translation','static'] if run==5 else ['rotation']):
        group=[r for r in rows if r['run']==run and r['condition']==condition];good={lane:[r for r in group if r['lane']==lane and r['result']['valid']] for lane in ['old','new']}
        common=set(r['case'] for r in good['old'])&set(r['case'] for r in good['new'])
        for lane in ['old','new']:
            valid=good[lane];paired=[r for r in valid if r['case'] in common];axes=np.array([r['result']['axis_camera'] for r in paired]);spread=None
            if len(axes)>1:
                mean=axes.mean(axis=0);spread=med([angle(a,mean) for a in axes]) if np.linalg.norm(mean)>1e-9 else None
            proj=[angle(r['result']['axis_camera'][:2],r['shaft']['direction'],True) for r in paired if r['shaft']['valid'] and np.linalg.norm(r['result']['axis_camera'][:2])>.2]
            summary.append(dict(run=run,condition=condition,lane=lane,planned=len(group)//2,valid=len(valid),paired=len(paired),
                rpm_median=med([r['result']['spin_rpm'] for r in valid]),reference_abs_deviation_paired_mean=float(np.mean([abs(r['result']['spin_rpm']-nominal) for r in paired])) if paired and run in [1,4] else None,
                shaft_projection_paired_median_deg=med(proj),shaft_count=len(proj),axis_spread_paired_median_deg=spread,
                reasons=dict(Counter(r['result']['failure_reason'] for r in group if r['lane']==lane and not r['result']['valid']))))
trial.old.save(OUT/'summary.json',summary)
print(json.dumps(summary,indent=2))
# Visual QA of fixed midpoint frames; missing geometry is not painted as a detection.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axes=plt.subplots(6,3,figsize=(12,20))
for route,ax in zip(routes,axes.flat):
    mid=route['end_us']-10000
    g=next((g for g in geometries if g['case']==route['case'] and g['mid_us']==mid),None)
    if g:
        raw=np.load(OUT/g['input_file'])['events'];c=np.array(g['anchor']);box=np.r_[np.floor(c-90).astype(int),np.ceil(c+90).astype(int)]
        ax.imshow(trial.p2.picture(raw,box),extent=[box[0],box[2],box[3],box[1]])
        ax.add_patch(plt.Circle(c,g['radius'],fill=False,color='gray',ls='--',lw=.7))
        pts=np.array(g['candidate']['points']).reshape(-1,2)
        if len(pts):ax.scatter(pts[:,0],pts[:,1],s=8,color='red')
        if g['fit']['valid']:ax.add_patch(plt.Circle(g['fit']['center'],g['radius'],fill=False,color='lime',lw=.8))
    else:ax.text(.5,.5,'No geometry anchor',ha='center',transform=ax.transAxes)
    ax.set_title(route['case']+' / '+route['new_route'],fontsize=9);ax.set_xticks([]);ax.set_yticks([])
fig.suptitle('Fixed new contexts | red actual contour pixels; green valid midpoint fit; dashed prior')
fig.tight_layout(rect=[0,0,1,.97]);fig.savefig(OUT/'geometry_review.png',dpi=110);plt.close(fig)
fmt=lambda v:'—' if v is None else f'{v:.2f}'
lines=['# 外侧轮廓固定参数扩展验证与回退','',
       '完成18个新增目标片段、36条路线状态；提取135个5ms几何快门。目标旋转区间与上一轮外侧前端数据的实际快门范围无重叠，仍来自同一组录制，非独立采集盲测。', '',
       '**结论：几何回退机制已验证，但外侧前端未显示稳定、全面的优势，保留可选实验方案，不升级生产默认。**', '',
       '## 主要结果', '',
       '- 新旧两路线均8/12旋转窗口可用，4个窗口没有几何支撑。新路线没有损失旧路线的有效输出，也没有恢复额外窗口。',
       '- 正反转共同六窗相对同档位参考的平均绝对偏差约13.50→13.33 RPM，差异很小；正转10.51→13.58变差，反转16.49→13.08改善，不能称稳定提精度。',
       '- 固定轴离散中位：run1为5.68→2.71度，run4为6.12→3.92度；run3仅两窗约1.03→1.06度，且有一窗回退，证据有限。',
       '- 纯平移0/4接受非零，两条路线均无新增误报；静止两窗均事件不足，未认证零速。',
       '- 杆投影run1只1个可评片段，run3也只1个且新路线回退旧结果，不能用该结果声称新前端斜轴改善已复现；run4三个共同支持窗中位3.41→3.75度。', '',
       '## 固定回退规则与实际使用', '',
       '在CMAX求解前：outer20完整几何覆盖可用则使用outer20；否则退到old20；old20不覆盖则old10；均不可用则拒绝。只检查几何时间覆盖与轮廓拟合有效性，不根据RPM/对比度/旋转求解是否成功回退。',
       '实际18片段：11个outer20、3个old20回退、4个unavailable。3个回退包含1个旋转与2个静止。相同路线显式复用同一个结果，不当作两次独立实验。本轮未实际触发old10分支，该分支有测试及此前实录证据，不能宣称本轮新增实录验证。', '',
       '|录制/条件|路线|有效/预定|RPM中位|共同参考偏差RPM|共同轴离散°|杆投影差°（支持数）|',
       '|---|---|---:|---:|---:|---:|---|']
for r in summary:lines.append(f"|{r['run']}/{r['condition']}|{r['lane']}|{r['valid']}/{r['planned']}|{fmt(r['rpm_median'])}|{fmt(r['reference_abs_deviation_paired_mean'])}|{fmt(r['axis_spread_paired_median_deg'])}|{fmt(r['shaft_projection_paired_median_deg'])}（{r['shaft_count']}）|")
lines += ['', '## 限制与后续', '',
          '旧新路线各自使用自己的中心和0.8r事件集，本次比较完整前端组合，不把差异归因为单一因素。样本小，旧中心先验和同录制片段有相关性。半径/杆方向仍是已知实验条件，固定相机/装置研究不能直接泛化飞行球。',
          '另拍同档位约246.91RPM只能提供弱参考；run3该方向未单独确认相同转速，未计算其参考偏差。二维杆投影与轴离散不是真实三维轴精度。轮廓与中值平滑使用未来几何支点，未验证在线时延。',
          '建议停止在这组小片段上反复选择参数，保留现有两方案和明确回退，用新录制或同步外部球心参考检验实际中心误差。若只做离线使用，可在实验目录显式选择新路线，不宣称已经可靠替代旧前端。', '',
          '## 核验与复现', '',
          '6项回退边界测试通过；135份几何输入、25份独立旋转输入的原始事件ID/时间范围/ROI/中心与抽样重建通过；18个路由由独立逐间隔检查验证；冻结源文件和求解前路由哈希未变。',
          '使用g1_spin/.venv/bin/python运行test_policy.py（pytest）、run.py、verify_report.py。run.py拒绝覆盖已有输出。',
          'policy.py定义回退，routes.json在求解前冻结；geometry.json、results.json和summary.json分别保存几何、逐条结果和汇总；verification.json保存核验；geometry_review.png为固定中点原始事件复核。']
# Derive the unique input count rather than relying on a written estimate.
lines=[line.replace('25份独立旋转输入',f'{len(checks)}份独立旋转输入') for line in lines]
(OUT/'report.md').write_text('\n'.join(lines)+'\n')
print('Report and geometry review saved')
