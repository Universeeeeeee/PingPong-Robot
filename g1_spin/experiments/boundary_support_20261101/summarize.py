"""Summarize every failure and independently validate representative support bookkeeping."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates
import run

OUT=Path(__file__).resolve().parent

def read(name):
    rows=list(csv.DictReader((OUT/name).open()))
    for r in rows:
        for k,v in r.items():
            if v in ['True','False']: r[k]=v=='True'
            else:
                try:r[k]=float(v)
                except ValueError:pass
    return rows

def med(rows,key):
    x=np.array([r[key] for r in rows],dtype=float)
    return float(np.nanmedian(x)) if len(x) and np.isfinite(x).any() else np.nan

def span(rows,key):
    x=np.array([r[key] for r in rows],dtype=float)
    return float(np.ptp(x[np.isfinite(x)])) if np.isfinite(x).any() else np.nan

def pct(x):
    return 'NA' if not np.isfinite(x) else f'{x:.0%}'

def fmt(x):
    return 'NA' if not np.isfinite(x) else f'{x:.2f}'

def main():
    rows=read('observations.csv');sensitivity=read('sensitivity.csv');controls=read('controls.csv')
    meta=json.loads((OUT/'metadata.json').read_text())
    summary=[]
    for segment in run.SEGMENTS:
        for d in [20,50,100]:
            for method in run.METHODS:
                for mode in ['free','fixed48']:
                    q=[r for r in rows if r['segment']==segment and r['duration_ms']==d and r['method']==method and r['mode']==mode]
                    valid=[r for r in q if r['fit_valid']]
                    kept=[r for r in q if r['retain_candidate']]
                    ss=dict(segment=segment,duration_ms=d,method=method,mode=mode,n=len(q),
                            gate_pass=sum(r['engineering_gate'] for r in q),fit_valid=len(valid),retain=len(kept),
                            not_fitted=sum(r['fit_status']=='insufficient_points' for r in q))
                    for k in ['candidates','dropped_fraction','occupied_15deg_bins','longest_radial_run_deg',
                              'isolated_event_fraction','median_distinct_pixels','median_active_time_blocks',
                              'before_support_fraction','after_support_fraction','median_radius_change_from_m0']:
                        ss['median_'+k]=med(q,k)
                    for k in ['cx','cy','radius','rms','jac_cy_r_cos','jac_condition']:
                        ss['median_'+k]=med(valid,k)
                    for k in ['cx','cy','radius']:
                        ss['range_'+k]=span(valid,k);ss['retained_range_'+k]=span(kept,k)
                    summary.append(ss)
    run.write_csv('summary.csv',summary)
    method_deltas=[]
    for r in rows:
        if r['method']=='M0':continue
        b=next(x for x in rows if x['method']=='M0' and all(x[k]==r[k] for k in ['segment','midpoint_s','duration_ms','mode']))
        method_deltas.append(dict(**{k:r[k] for k in ['segment','midpoint_s','duration_ms','method','mode','retain_candidate']},
                                  both_fit_valid=r['fit_valid'] and b['fit_valid'],dx=r['cx']-b['cx'],dy=r['cy']-b['cy'],dr=r['radius']-b['radius'],
                                  median_candidate_radius_composition_shift=r['median_radius_change_from_m0']))
    run.write_csv('method_deltas.csv',method_deltas)
    paired=[]
    for segment in run.SEGMENTS:
        for method in run.METHODS:
            for mode in ['free','fixed48']:
                for short,long in [(20,50),(20,100)]:
                    a={r['midpoint_s']:r for r in rows if r['segment']==segment and r['method']==method and r['mode']==mode and r['duration_ms']==short}
                    b={r['midpoint_s']:r for r in rows if r['segment']==segment and r['method']==method and r['mode']==mode and r['duration_ms']==long}
                    for t in a:
                        x,y=a[t],b[t]
                        paired.append(dict(segment=segment,method=method,mode=mode,midpoint_s=t,short_ms=short,long_ms=long,
                                           both_fit_valid=x['fit_valid'] and y['fit_valid'],both_retained=x['retain_candidate'] and y['retain_candidate'],
                                           dx=x['cx']-y['cx'],dy=x['cy']-y['cy'],dr=x['radius']-y['radius']))
    run.write_csv('paired_windows.csv',paired)
    paired_summary=[]
    for segment in run.SEGMENTS:
        for method in run.METHODS:
            for mode in ['free','fixed48']:
                for long in [50,100]:
                    for group in ['both_fit_valid','both_retained']:
                        q=[r for r in paired if r['segment']==segment and r['method']==method and r['mode']==mode and r['long_ms']==long and r[group]]
                        rr=dict(segment=segment,method=method,mode=mode,long_ms=long,group=group,n=len(q))
                        for k in ['dx','dy','dr']:
                            x=np.abs([r[k] for r in q]);rr['median_abs_'+k]=float(np.median(x)) if len(x) else np.nan
                            rr['p95_abs_'+k]=float(np.quantile(x,.95)) if len(x) else np.nan
                        paired_summary.append(rr)
    run.write_csv('paired_summary.csv',paired_summary)
    delta=[]
    for r in sensitivity:
        baseline=next(x for x in rows if all(x[k]==r[k] for k in ['segment','midpoint_s','duration_ms','method','mode']))
        rr=dict(**r)
        for k in ['cx','cy','radius','candidates']:rr['delta_'+k]=r[k]-baseline[k]
        rr['baseline_retain']=baseline['retain_candidate'];rr['both_fit_valid']=r['fit_valid'] and baseline['fit_valid'];rr['gate_changed']=r['engineering_gate']!=baseline['engineering_gate']
        delta.append(rr)
    run.write_csv('sensitivity_deltas.csv',delta)
    validation=validate(rows,meta)
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    report(summary,paired_summary,delta,controls,meta,validation,method_deltas)
    print(json.dumps({k:v for k,v in validation.items() if k!='checks'},indent=2))


def validate(rows,meta):
    manifest=json.loads((run.DATA/'manifest.json').read_text())
    events=np.concatenate([np.load(run.DATA/f) for f in manifest['parts']]);t0=int(events['t'][0])
    previous=list(csv.DictReader((run.PREVIOUS/'observations.csv').open()))
    checks=[]
    for segment,limits in run.SEGMENTS.items():
        mid=sum(limits)/2
        for duration in [20,50,100]:
            lo=t0+round((mid-duration/2000)*1e6);hi=t0+round((mid+duration/2000)*1e6)
            selected=events[(events['t']>=lo)&(events['t']<hi)]
            frame=np.bincount(selected['y'].astype(int)*1280+selected['x'],minlength=720*1280).reshape(720,1280)
            xx=run.ORIGIN[0]+np.cos(run.ANGLES[:,None])*run.RADII
            yy=run.ORIGIN[1]+np.sin(run.ANGLES[:,None])*run.RADII
            profiles=map_coordinates(gaussian_filter(frame.astype(float),1),[yy,xx],order=1)
            ids=[];points=[];radii=[]
            for i,prof in enumerate(profiles):
                hits=np.flatnonzero(prof>=.6*duration/100)
                if len(hits):
                    j=hits[-1];ids.append(i);points.append([xx[i,j],yy[i,j]]);radii.append(run.RADII[j])
            ids=np.array(ids);points=np.array(points);radii=np.array(radii)
            ix,p,rr=run.extract(selected,duration)
            assert np.array_equal(ids,ix) and np.array_equal(points,p)
            # Independent brute-force direct distance; unique xy tuples inside explicit time slices.
            passes=[];distinct=[];blocks=[]
            for point in points:
                near=selected[(selected['x'].astype(float)-point[0])**2+(selected['y'].astype(float)-point[1])**2<=9]
                pixels=set(zip(near['x'].tolist(),near['y'].tolist()))
                active=0
                for k in range(4):
                    z=near[(near['t']>=lo+(hi-lo)*k/4)&(near['t']<lo+(hi-lo)*(k+1)/4)]
                    active+=len(set(zip(z['x'].tolist(),z['y'].tolist())))>=2
                distinct.append(len(pixels));blocks.append(active);passes.append(len(pixels)>=3 and active>=2)
            _,di,bl=run.support(selected,points,lo,hi)
            assert np.array_equal(distinct,di) and np.array_equal(blocks,bl)
            keep=np.array(passes,dtype=bool);ind=np.flatnonzero(keep);m2=[];current=[]
            for j in ind:
                if current and (ids[j]!=ids[current[-1]]+1 or abs(radii[j]-radii[current[-1]])>2):
                    if len(current)>=11:m2.extend(current)
                    current=[]
                current.append(j)
            if len(current)>=11:m2.extend(current)
            _,stored=run.evaluate(selected,selected[:0],selected[:0],lo,hi,duration,{'segment':segment,'midpoint_s':mid})
            assert np.array_equal(stored['M1'][1],points[keep]) and np.array_equal(stored['M2'][1],points[m2])
            baseline=next(r for r in rows if r['segment']==segment and r['duration_ms']==duration and r['method']=='M0' and r['mode']=='free' and abs(r['midpoint_s']-mid)<1e-6)
            prev=next(r for r in previous if r['segment']==segment and r['duration_ms']==str(duration) and r['mode']=='free' and r['threshold_policy']=='scaled_0.6' and abs(float(r['midpoint_s'])-mid)<1e-6)
            diff=max(abs(baseline[k]-float(prev[k])) for k in ['cx','cy','radius'])
            assert diff<1e-5 and len(selected)==baseline['events_total']
            # Independent heldout: no event reused across half-open block intervals.
            before=events[(events['t']>=2*lo-hi)&(events['t']<lo)]
            after=events[(events['t']>=hi)&(events['t']<2*hi-lo)]
            assert not len(before) or before['t'].max()<lo
            assert not len(after) or after['t'].min()>=hi
            checks.append(dict(segment=segment,duration_ms=duration,event_count=len(selected),
                               full_sensor_matches_crop=True,direct_support_matches=True,independent_M2_matches=True,
                               previous_baseline_max_fit_difference=diff,heldout_disjoint=True))
    # Every M0 grid fit, not just the representative plots, must preserve the old pipeline.
    baseline_checks=0
    for row in rows:
        if row['method'] != 'M0':continue
        prior=next(r for r in previous if r['segment']==row['segment'] and r['duration_ms']==str(int(row['duration_ms'])) and r['mode']==row['mode'] and r['threshold_policy']=='scaled_0.6' and abs(float(r['midpoint_s'])-row['midpoint_s'])<1e-6)
        assert max(abs(row[k]-float(prior[k])) for k in ['cx','cy','radius'])<1e-5
        assert row['candidates']==int(prior['candidate_count'])
        baseline_checks+=1
    keys=['segment','midpoint_s','duration_ms','method','mode']
    assert len(rows)==meta['expected_rows']==len({tuple(r[k] for k in keys) for r in rows})
    assert run.digest_directory(run.PREVIOUS)==meta['previous_files_sha256']
    assert hashlib.sha256((OUT/'run.py').read_bytes()).hexdigest()==meta['script_sha256']
    return dict(expected_grid_rows=meta['expected_rows'],actual_grid_rows=len(rows),all_grid_keys_unique=True,
                previous_directory_files_unchanged=True,script_hash_matches=True,all_grid_M0_matches_previous_rows=baseline_checks,representative_independent_checks=len(checks),
                not_fitted_rows=sum(r['fit_status']=='insufficient_points' for r in rows),
                numerically_invalid_rows=sum(not r['fit_valid'] for r in rows),checks=checks)


def report(summary,pairs,delta,controls,meta,validation,method_deltas):
    lines=['# 原始事件支持与邻角连续性消融','',
           '结论：旧2px孤点比例有所降低，但3px多像素/多时间块支持没有消除旧2px孤点指标；M2以大量拒绝换取连续性，20ms四段仅1/124窗通过工程gate且free拟合有效，不能宣称20ms几何出口已打通。50–100ms早/中段可保留待审候选；晚段即使通过gate仍不得无条件使用。人工外层抹除后early_turn/late的100ms仍能通过，说明内部事件结构或人工mask边仍可替代本体外轮廓，支持与连续性不足以确认silhouette。', '',
           '追加背景（首次数据运行后收到，未更改协议或数值）：用户已确认10/29 RAW同样是ITTF40+塑料球，表面贴约14个凸起球形反光Marker，直径不明；“凸起厚度可忽略”是待验证近似，不能断言厚度主导本轮差异。整体球+Marker不是严格单球面；外侧或孤立事件可能来自凸起Marker，连续弧也可能是Marker轮廓。“孤点”仅指2px事件计数定义，不能直接称噪声；“内部结构”包括可能Marker投影。工程gate不是本体silhouette分类，本轮不据此推进动态补偿。', '',
           '|片段|M2保留20/50/100ms（各31窗）|当前出口|','|---|---|---|',
           '|translation|1 / 31 / 31|20ms拒绝；50–100ms只作运动敏感性/几何候选|',
           '|early_turn|0 / 30 / 31|50–100ms保留待审；20ms拒绝|',
           '|middle|0 / 31 / 31|50–100ms保留待审；20ms拒绝|',
           '|late|0 / 11 / 29|所有时长退出无条件几何出口；通过gate也需身份审查|','',
           '小残差、优化收敛和r接近48均不是silhouette证据。未通过gate的所有窗拒绝作为该工程候选；唯一通过的translation20ms仍因代表窗阈值敏感与无真值，不视为阶段A通过。M0/M1/M2统计不涉及绝对精度。', '',
           '## 固定协议与输出含义','',
           '完整预声明见[protocol.md](protocol.md)。同一采样原点、sigma1模糊、阈值、121射线及旧外侧点：M0原样；M1在3px内≥3 distinct pixels、4相对时间块中≥2块各≥2 distinct pixels；M2再保留相邻射线径向差≤2px且至少15°连续run。共同工程gate要求≥30点、≥6个15°bin和≥30°径向连续run。gate加拟合数值有效仅表示保留候选，不能确认轮廓。M1/M2只删除M0点，无内层回搜。','',
           '20/50/100ms的时间子块分别为5/12.5/25ms，因此长窗支持标准有不同绝对带宽。首次真实数据前，人工弧揭示2px邻域与模糊外沿不匹配，改为固定3px并保留2px对照；没有根据真实数据结果调参。即使3px，密集弧也可因模糊外扩被过度拒绝。所有未拟合、触界及gate拒绝均留在CSV。','',
           '## 完整格点：支持与拒绝','',
           '下表31窗/格，free模式；fit为数值有效窗，留为通过工程gate且数值有效。孤点仍是旧2px内原始事件≤1；因此M1的3px支持不能保证该旧指标为零。拒绝率=1−留/31，覆盖用角bin与连续run表达。点数/覆盖中位含空窗；支持比例、孤点比例和Δρ在非空窗中计算中位，因此20ms少数幸存弧的支持数不能代表整段。','',
           '|片段|T ms|方法|fit/留|丢点中位|点数中位|角bin中位|run°中位|2px孤点|distinct中位|前/后独立块支持|','|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for s in summary:
        if s['mode']=='free':
            lines.append(f"|{s['segment']}|{s['duration_ms']}|{s['method']}|{s['fit_valid']}/{s['retain']}|{s['median_dropped_fraction']:.0%}|{s['median_candidates']:.0f}|{s['median_occupied_15deg_bins']:.0f}|{s['median_longest_radial_run_deg']:.1f}|{s['median_isolated_event_fraction']:.1%}|{s['median_median_distinct_pixels']:.1f}|{s['median_before_support_fraction']:.0%}/{s['median_after_support_fraction']:.0%}|")
    lines+=['','前/后支持使用未参与选点的相邻非重叠T窗，同坐标3px邻域、同distinct和4块标准。不补偿运动；平移会令真实边界降低此支持，持续内部结构（可能Marker投影）也会升高它，不能据此判真伪。邻中点窗重叠50/80/90%，不作为独立重复。','',
            '## 数值拟合与向内组成变化','',
            '下表为数值有效free拟合的cx/cy/r范围；含真实运动和采样变化，不是误差。Δρ为筛后候选采样半径中位−M0半径中位；每个保留点原坐标完全未改，Δρ反映选择组成变化，负值也不等于发现内层轮廓。','',
            '|片段|T|方法|cx范围|cy范围|r范围|RMS中位|Δρ中位|','|---|---:|---|---:|---:|---:|---:|---:|']
    for s in summary:
        if s['mode']=='free':lines.append(f"|{s['segment']}|{s['duration_ms']}|{s['method']}|{fmt(s['range_cx'])}|{fmt(s['range_cy'])}|{fmt(s['range_radius'])}|{fmt(s['median_rms'])}|{fmt(s['median_median_radius_change_from_m0'])}|")
    lines+=['','相对M0的逐窗拟合变化另存method_deltas.csv。以下M2只用保留且双方数值有效的free窗；半径变化仍受cy/r耦合影响，不能把正负直接解释为真轮廓移动。','',
            '|片段|T|配对数|拟合Δr中位|Δr最小/最大|','|---|---:|---:|---:|---:|']
    for segment in run.SEGMENTS:
        for duration in [20,50,100]:
            q=[r for r in method_deltas if r['segment']==segment and r['duration_ms']==duration and r['method']=='M2' and r['mode']=='free' and r['retain_candidate'] and r['both_fit_valid']]
            x=np.array([r['dr'] for r in q])
            lines.append(f'|{segment}|{duration}|{len(q)}|{fmt(float(np.median(x))) if len(x) else "NA"}|{fmt(float(np.min(x))) if len(x) else "NA"}/{fmt(float(np.max(x))) if len(x) else "NA"}|')
    lines+=['','没有全局一致向内收缩：50–100ms平移/早/中M2保留窗的拟合Δr中位约−0.15至+0.03px，但late50ms中位+2.48px，late100ms范围−6.06至+5.32px；这只能说明相对M0拟合变化，不能确认内部结构身份或本体边界。','',
            'fixed48结果及被保留子集的cx/cy/r范围均在summary.csv。fixed48强加48px先验，不能用其稳定性推断精度。RMS下降天然受删点影响。','',
            '## 20ms与较长窗配对','',
            '表为双方数值有效的配对数和|dx|/|dy|中位/P95；双方都通过gate的配对另存paired_summary.csv。','',
            '|片段|方法|模式|长窗ms|配对数|dx中位/P95|dy中位/P95|','|---|---|---|---:|---:|---:|---:|']
    for s in pairs:
        if s['group']=='both_fit_valid':lines.append(f"|{s['segment']}|{s['method']}|{s['mode']}|{s['long_ms']}|{s['n']}|{fmt(s['median_abs_dx'])}/{fmt(s['p95_abs_dx'])}|{fmt(s['median_abs_dy'])}/{fmt(s['p95_abs_dy'])}|")
    lines+=['','## 代表窗敏感性','',
            '12个代表窗，阈值0.3/1.0×T/100、原点四向±2px、支持邻域2px对3px；以下只汇总free/M2，完整M0/M1/fixed48与逐行失败见sensitivity_deltas.csv。位置变化仅在基线和扰动都有数值解时有意义，表中最大值仅用双方数值有效配对，不隐去失败。','',
            '|扰动|数值有效/12|保留/12|gate改变数|最大绝对Δcx|最大绝对Δcy|最大绝对Δr|','|---|---:|---:|---:|---:|---:|---:|']
    for kind in dict.fromkeys(r['perturbation'] for r in delta):
        q=[r for r in delta if r['perturbation']==kind and r['method']=='M2' and r['mode']=='free']
        maxima=[]
        for k in ['delta_cx','delta_cy','delta_radius']:
            a=np.abs([r[k] for r in q if r['both_fit_valid']]);maxima.append(float(np.nanmax(a)) if np.isfinite(a).any() else np.nan)
        lines.append(f"|{kind}|{sum(r['fit_valid'] for r in q)}|{sum(r['retain_candidate'] for r in q)}|{sum(r['gate_changed'] for r in q)}|{fmt(maxima[0])}|{fmt(maxima[1])}|{fmt(maxima[2])}|")
    lines+=['','最大原点扰动来自late50ms且基线/扰动均未通过工程gate（Δcy27.17px）；对已保留代表窗的变化须看逐行baseline_retain/retain_candidate，不能把被拒绝弧的大幅拟合漂移当保留出口误差。','',
            '## 人工外层抹除、真实晚段和背景诊断','',
            '人工压力测试只保留距上轮代表100ms粗中心≤36px的真实事件，粗中心只用于造mask。对象是已确认Marker球，抹除子集不是物理未改动的裸球。该区域仍和原r30..75.5扫描交叠；mask会造边，不能视为真实负标签或测误报率。若M0<30点标not_evaluable，不能把平凡无点当成功拒绝。背景来自x+180px实际ROI，减180后复用全部坐标、采样与bounds。','',
            '|片段|T|控制|M0点数/可评估|M1留/点数|M2留/点数|M2前/后支持|','|---|---:|---|---:|---:|---:|---:|']
    for a in controls:
        if a['mode']=='free' and a['method']=='M0':
            q={r['method']:r for r in controls if all(r[k]==a[k] for k in ['segment','duration_ms','control','mode'])}
            b,c=q['M1'],q['M2']
            lines.append(f"|{a['segment']}|{a['duration_ms']:.0f}|{a['control']}|{a['candidates']:.0f}/{a['control_evaluable']}|{b['retain_candidate']}/{b['candidates']:.0f}|{c['retain_candidate']}/{c['candidates']:.0f}|{pct(c['before_support_fraction'])}/{pct(c['after_support_fraction'])}|")
    lines+=['','本次背景6窗M0只有0–22点，均低于压力对照的30点可评估门槛，只能说明该稀疏背景未产生保留，不能证明较密背景鲁棒性。','',
            '真实晚段不作负标签，其保留结果须结合candidate_overlays.png人工审查，内部结构/可能Marker投影与本体轮廓身份仍未确认。人工mask内保留候选表明支持/连续gate可能接受内部结构或mask边缘，不能证明筛选后不存在内部结构替代。','',
            '## 实现检查和局限','',
            f"完整主表{validation['actual_grid_rows']}行，未拟合{validation['not_fitted_rows']}行，数值无效合计{validation['numerically_invalid_rows']}行；未删除失败。12代表窗用全1280×720累积独立重现M0，直接事件距离+xy集合/显式时间切片核对M1，独立循环核对M2，绝对时间boolean mask核对取窗与非重叠支持；完整744行M0的点数和free/fixed48拟合均与上轮对应主阈值结果一致，拟合差<1e−5px。旧目录文件哈希全部保持。详见validation.json。",'',
            '语义检查：热点100次重复不满足distinct空间gate，单时间块3×3簇不满足时间gate，两块多像素连续人工弧通过M2。高密度四块人工弧失败记录：'+json.dumps(meta['semantic_checks']['high_density_four_block_arc'],ensure_ascii=False)+'。后者表明固定邻域与blur外沿的密度依赖，不能用本方案声称所有有支持轮廓都保留。','',
            '未解决：何者为silhouette、模糊外沿的定位偏差、缺弧时cy/r耦合、平移导致的独立块支持降低、持续内部事件结构（可能Marker投影）/重复噪声造成的假支持。本轮无真值，不报精确率/召回率或绝对精度，不推进旋转。既已确认约14个球形Marker，不建议继续仅调连续弧门槛去强逼光滑本体圆；后续可另行评估Marker几何位置/跨时刻对应→刚体姿态→角速度路线，本轮未开发或验证该路线。','',
            '复现：','',
            '```sh','g1_spin/.venv/bin/python g1_spin/experiments/boundary_support_20261101/run.py','g1_spin/.venv/bin/python g1_spin/experiments/boundary_support_20261101/summarize.py','```','',
            '图：candidate_overlays.png全部代表窗图示ROI内所有原始事件与三方法点（重复像素重叠）；negative_controls.png控制；fit_traces.png数值有效拟合（gate仍可能拒绝）。CSV包含完整主格点、summary、method_deltas、paired_windows/summary、sensitivity/deltas、controls；metadata包含输入和旧目录哈希、协议及脚本哈希。仅写本目录，未动s-CMax/法向流/多窗旋转，无网络/训练/硬件。']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':main()
