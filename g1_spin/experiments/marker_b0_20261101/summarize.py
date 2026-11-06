"""B0 per-seed outcomes, model diagnostics, cross-window checks and independent audits."""
import csv
import hashlib
import json
from pathlib import Path
import collections

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter

import run
OUT=Path(__file__).resolve().parent
KINDS={'anchor':('anchor_x','anchor_y'), 'combined_centroid':('combined_centroid_x','combined_centroid_y'),
       'on_centroid':('on_centroid_x','on_centroid_y'),'off_centroid':('off_centroid_x','off_centroid_y'),
       'circle':('circle_cx','circle_cy'),'ellipse':('ellipse_cx','ellipse_cy')}


def read(name):
    rows=list(csv.DictReader((OUT/name).open()))
    for r in rows:
        for k,v in r.items():
            if v in ['True','False']:r[k]=v=='True'
            elif not v:r[k]=np.nan
            else:
                try:r[k]=float(v)
                except ValueError:pass
    return rows


def point(r,kind):
    if not r['position_valid']:return None
    x,y=KINDS[kind]
    v=np.array([r.get(x,np.nan),r.get(y,np.nan)])
    return v if np.isfinite(v).all() else None


def stats(v):
    v=np.array(v,dtype=float);v=v[np.isfinite(v)]
    return dict(n=len(v),median=float(np.median(v)) if len(v) else np.nan,
                p95=float(np.quantile(v,.95)) if len(v) else np.nan,max=float(np.max(v)) if len(v) else np.nan)


def f(x):return 'NA' if not np.isfinite(x) else f'{x:.2f}'


def main():
    rows=read('observations.csv');summary=[];cross=[];comparisons=[]
    for segment in run.raw.SEGMENTS:
        for feature in run.SEEDS[segment]:
            for duration in [2,5]:
                q=[r for r in rows if r['segment']==segment and r['feature']==feature and r['duration_ms']==duration]
                valid=[r for r in q if r['position_valid']]
                refusals=[r for r in q if str(r['status']).startswith('rejected_')]
                summary.append(dict(segment=segment,feature=feature,duration_ms=duration,frames=len(q),valid_with_manual_seed=len(valid),
                                    automatic_accepted=sum(r['status']=='accepted' for r in q),not_attempted=sum(r['status']=='not_attempted_after_loss' for r in q),
                                    valid_start_s=min(r['start_s'] for r in valid),valid_end_s=max(r['end_s'] for r in valid),
                                    seed_block_midpoint_s=next(r['midpoint_s'] for r in q if r['status']=='manual_seed'),
                                    backward_termination=';'.join(f"frame{int(r['frame'])}:{r['status']}" for r in refusals if r['frame']<len(q)//2) or 'segment_boundary',
                                    forward_termination=';'.join(f"frame{int(r['frame'])}:{r['status']}" for r in refusals if r['frame']>len(q)//2) or 'segment_boundary',
                                    circle_valid=sum(r.get('circle_valid',False) is True for r in valid),
                                    ellipse_valid=sum(r.get('ellipse_valid',False) is True for r in valid),
                                    min_local_distinct_pixels=min(r['local_distinct_pixels'] for r in valid),
                                    median_local_distinct_pixels=float(np.median([r['local_distinct_pixels'] for r in valid]))))
                for r in valid:
                    for a,b in [('on_centroid','off_centroid'),('combined_centroid','anchor'),('circle','combined_centroid'),('ellipse','combined_centroid'),('circle','ellipse')]:
                        pa,pb=point(r,a),point(r,b)
                        if pa is not None and pb is not None:
                            delta=pa-pb
                            comparisons.append(dict(segment=segment,feature=feature,duration_ms=duration,frame=r['frame'],midpoint_s=r['midpoint_s'],
                                                    comparison=a+'_minus_'+b,dx=delta[0],dy=delta[1],distance=float(np.linalg.norm(delta))))
            a=[r for r in rows if r['segment']==segment and r['feature']==feature and r['duration_ms']==2]
            b=[r for r in rows if r['segment']==segment and r['feature']==feature and r['duration_ms']==5]
            times=np.array([r['midpoint_s'] for r in a])
            for kind in KINDS:
                aligned=[]
                for r in b:
                    p5=point(r,kind)
                    if p5 is None:continue
                    t=r['midpoint_s'];right=int(np.searchsorted(times,t));left=right-1
                    if left<0 or right>=len(a):continue
                    p0,p1=point(a[left],kind),point(a[right],kind)
                    if p0 is None or p1 is None:continue
                    # Only adjacent original 2ms indices, never bridge a rejected/lost block.
                    assert a[right]['frame']-a[left]['frame']==1
                    w=(t-times[left])/(times[right]-times[left]);p2=p0*(1-w)+p1*w
                    delta=p2-p5
                    aligned.append(dict(segment=segment,feature=feature,position_kind=kind,midpoint_s=t,
                                        two_ms_left_frame=a[left]['frame'],two_ms_right_frame=a[right]['frame'],interpolation_weight=w,
                                        x_2ms_interpolated=p2[0],y_2ms_interpolated=p2[1],x_5ms=p5[0],y_5ms=p5[1],
                                        raw_dx=delta[0],raw_dy=delta[1],raw_distance=float(np.linalg.norm(delta))))
                if aligned:
                    nominal=sum(run.raw.SEGMENTS[segment])/2
                    reference=min(aligned,key=lambda r:abs(r['midpoint_s']-nominal))
                    offset=np.array([reference['raw_dx'],reference['raw_dy']])
                    for rr in aligned:
                        movement=np.array([rr['raw_dx'],rr['raw_dy']])-offset
                        rr.update(common_reference_s=reference['midpoint_s'],fixed_initial_offset_x=offset[0],fixed_initial_offset_y=offset[1],
                                  relative_displacement_dx=movement[0],relative_displacement_dy=movement[1],relative_displacement_distance=float(np.linalg.norm(movement)))
                    cross.extend(aligned)
    run.write_csv('track_summary.csv',summary);run.write_csv('position_comparisons.csv',comparisons)
    run.write_csv('cross_window_positions.csv',cross)
    distributions=[]
    for segment in run.raw.SEGMENTS:
        for feature in run.SEEDS[segment]:
            for duration in [2,5]:
                q=[r for r in comparisons if r['segment']==segment and r['feature']==feature and r['duration_ms']==duration]
                for name in dict.fromkeys(r['comparison'] for r in q):
                    z=[r for r in q if r['comparison']==name]
                    offsets=np.array([[r['dx'],r['dy']] for r in z])
                    diameter=float(np.max(np.linalg.norm(offsets[:,None]-offsets[None,:],axis=2))) if len(z) else np.nan
                    distributions.append(dict(segment=segment,feature=feature,duration_ms=duration,diagnostic=name,
                                              **stats([r['distance'] for r in z]),offset_change_diameter_px=diameter))
            for kind in KINDS:
                q=[r for r in cross if r['segment']==segment and r['feature']==feature and r['position_kind']==kind]
                for measure in ['raw_distance','relative_displacement_distance']:
                    distributions.append(dict(segment=segment,feature=feature,duration_ms='2_vs_5',diagnostic=kind+'_'+measure,
                                              **stats([r[measure] for r in q])))
    run.write_csv('diagnostic_distributions.csv',distributions)
    overlap=[]
    for segment in run.raw.SEGMENTS:
        for duration in [2,5]:
            for frame in range(100//duration):
                q=[r for r in rows if r['segment']==segment and r['duration_ms']==duration and r['frame']==frame and r['position_valid']]
                for i,a in enumerate(q):
                    for b in q[i+1:]:
                        delta=np.abs(point(a,'anchor')-point(b,'anchor'));intersection=float(np.prod(np.maximum(0,27-delta)))
                        overlap.append(dict(segment=segment,duration_ms=duration,frame=frame,midpoint_s=a['midpoint_s'],
                                            feature_a=a['feature'],feature_b=b['feature'],anchor_distance_px=float(np.linalg.norm(delta)),
                                            patch_iou=intersection/(2*27*27-intersection)))
    run.write_csv('cross_feature_overlap.csv',overlap)
    meta=json.loads((OUT/'metadata.json').read_text())
    meta['matching_image_roi_inclusive']=[run.X0,run.X0+run.W-1,run.Y0,run.Y0+run.H-1]
    meta['future_fields_meaning']='Next non-overlapping block in PROPAGATION direction; direction=-1 means chronologically earlier.'
    meta['diagnostic_additions_only']='Cross-feature patch overlap and actual-frame review added after fixed tracking run; no gate or observations changed.'
    (OUT/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    validation=validate(rows)
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    plot_cross(cross);plot_failure_cases(rows);plot_actual_frames(rows);report(rows,summary,distributions,validation,overlap)
    print(json.dumps({k:v for k,v in validation.items() if k!='independent_score_checks'},indent=2))


def validate(rows):
    meta=json.loads((OUT/'metadata.json').read_text());events,_=run.raw.load_events();times=events['t'].astype(np.int64);t0=int(times[0])
    assert len(rows)==630==len({(r['segment'],r['duration_ms'],r['feature'],r['frame']) for r in rows})
    assert hashlib.sha256((OUT/'run.py').read_bytes()).hexdigest()==meta['script_sha256']
    assert hashlib.sha256((OUT/'protocol.md').read_bytes()).hexdigest()==meta['protocol_sha256']
    score_checks=[];coverage_checks=0
    for segment,(start,end) in run.raw.SEGMENTS.items():
        for duration in [2,5]:
            for feature in run.SEEDS[segment]:
                q=sorted([r for r in rows if r['segment']==segment and r['duration_ms']==duration and r['feature']==feature],key=lambda r:r['frame'])
                reference=len(q)//2
                for direction in [-1,1]:
                    lost=False
                    for k in range(reference+direction,len(q) if direction>0 else -1,direction):
                        r=q[k]
                        if lost:assert r['status']=='not_attempted_after_loss' and not r['position_valid']
                        if str(r['status']).startswith('rejected_'):lost=True
                # One accepted and one rejected pair per path where available: direct patch cosine.
                for want in ['accepted','rejected']:
                    candidates=[r for r in q if r['status']=='accepted' if want=='accepted'] if want=='accepted' else [r for r in q if str(r['status']).startswith('rejected_')]
                    if not candidates:continue
                    r=min(candidates,key=lambda x:abs(x['frame']-reference));prev=q[int(r['source_frame'])]
                    def selected(frame):
                        lo=t0+round(start*1e6)+int(frame)*duration*1000;hi=lo+duration*1000
                        return events[(times>=lo)&(times<hi)]
                    source=selected(prev['frame']);target=selected(r['frame'])
                    # Independently accumulate unique pixels via explicit loops, then smooth.
                    def independent_image(e):
                        im=np.zeros((2,run.H,run.W))
                        for x,y,p in set(zip(e['x'].tolist(),e['y'].tolist(),(e['p']>0).tolist())):
                            if run.X0<=x<run.X0+run.W and run.Y0<=y<run.Y0+run.H:im[int(p),y-run.Y0,x-run.X0]=1.
                        return np.stack([gaussian_filter(channel,.8) for channel in im])
                    source_image,target_image=independent_image(source),independent_image(target)
                    source_center=np.array([prev['anchor_x'],prev['anchor_y']]);target_center=np.array([r['proposed_x'],r['proposed_y']])
                    def direct_patch(im,center):
                        x,y=np.rint(center).astype(int)-[run.X0,run.Y0];return im[:,y-13:y+14,x-13:x+14]
                    a=direct_patch(source_image,source_center);b=direct_patch(target_image,target_center)
                    direct_score=float(np.sum(a*b)/(np.linalg.norm(a)*np.linalg.norm(b)))
                    assert abs(direct_score-r['score'])<1e-10
                    # Direct scoring every displacement also checks best sign and second-peak exclusion.
                    score_map=np.zeros((25,25))
                    for iy,dy in enumerate(range(-12,13)):
                        for ix,dx in enumerate(range(-12,13)):
                            p=direct_patch(target_image,source_center+[dx,dy]);den=np.linalg.norm(a)*np.linalg.norm(p)
                            score_map[iy,ix]=np.sum(a*p)/den if den>0 else 0
                    iy,ix=np.unravel_index(np.argmax(score_map),score_map.shape)
                    assert np.array_equal(source_center+[ix-12,iy-12],target_center)
                    score_map[max(0,iy-3):iy+4,max(0,ix-3):ix+4]=-np.inf
                    sy,sx=np.unravel_index(np.argmax(score_map),score_map.shape)
                    assert abs(score_map[sy,sx]-r['second_score'])<1e-10
                    score_checks.append(dict(segment=segment,duration_ms=duration,feature=feature,status=r['status'],frame=r['frame'],
                                             absolute_time_mask=True,direct_cosine_matches=True,all_displacements_and_second_peak_match=True))
                for r in q:
                    if not r['position_valid']:continue
                    for kind in ['circle','ellipse']:
                        if r.get(kind+'_valid',False) is True:
                            assert r[kind+'_occupied_30deg_bins']>=8 and r[kind+'_max_empty_angle_deg']<=90
                            assert r[kind+'_median_abs_residual']<=1.5 and r[kind+'_rms']<=2.5
                            coverage_checks+=1
    return dict(expected_rows=630,actual_rows=len(rows),all_unique=True,run_and_protocol_hash_match=True,
                nonreconnection_verified_all_36_directions=True,independently_rescored_pairs=len(score_checks),
                valid_geometry_gate_rows_checked=coverage_checks,semantic_checks=meta['semantic_checks'],
                status_counts=dict(collections.Counter(r['status'] for r in rows)),independent_score_checks=score_checks,
                summarize_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


def plot_cross(cross):
    fig,axes=plt.subplots(3,3,figsize=(14,9))
    for i,segment in enumerate(run.raw.SEGMENTS):
        for j,feature in enumerate(run.SEEDS[segment]):
            ax=axes[i,j];q=[r for r in cross if r['segment']==segment and r['feature']==feature and r['position_kind']=='anchor']
            ax.plot([r['midpoint_s'] for r in q],[r['raw_distance'] for r in q],'o-',ms=3,label='raw anchor difference')
            ax.plot([r['midpoint_s'] for r in q],[r['relative_displacement_distance'] for r in q],'s-',ms=3,label='after ONE fixed initial offset')
            ax.set_title(f'{segment} {feature}; paired n={len(q)}');ax.set_ylabel('2ms vs 5ms distance (px)');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=7)
    for ax in axes[-1]:ax.set_xlabel('Aligned time (s); no interpolation across loss')
    fig.suptitle('Cross-window diagnostic; differing seed times and cumulative patch drift are not ground truth',fontsize=12)
    fig.tight_layout(rect=(0,0,1,.96));fig.savefig(OUT/'cross_window_differences.png',dpi=180);plt.close(fig)


def plot_failure_cases(rows):
    events,_=run.raw.load_events();times=events['t'].astype(np.int64);t0=int(times[0]);fig,axes=plt.subplots(6,3,figsize=(12,19))
    for i,(segment,duration) in enumerate((s,d) for s in run.raw.SEGMENTS for d in [2,5]):
        failures=[r for r in rows if r['segment']==segment and r['duration_ms']==duration and str(r['status']).startswith('rejected_')]
        if not failures:
            for ax in axes[i]:ax.text(.5,.5,f'{segment} {duration}ms: no rejection\nAll three paths reached segment boundaries',ha='center',va='center',transform=ax.transAxes);ax.axis('off')
            continue
        r=min(failures,key=lambda r:(abs(r['frame']-r['reference_frame']),r['feature'],r['frame']))
        source=next(x for x in rows if x['segment']==segment and x['duration_ms']==duration and x['feature']==r['feature'] and x['frame']==r['source_frame'])
        start=run.raw.SEGMENTS[segment][0]
        source_center=np.array([source['anchor_x'],source['anchor_y']]);best=np.array([r['proposed_x'],r['proposed_y']]);second=np.array([r['second_x'],r['second_y']])
        for ax,frame,center,label in zip(axes[i],[source['frame'],r['frame'],r['frame']],[source_center,best,second],['source accepted patch','rejected best proposal','alternative peak proposal']):
            lo=t0+round(start*1e6)+int(frame)*duration*1000;hi=lo+duration*1000
            e=events[(times>=lo)&(times<hi)];run.raw.draw_raw(ax,run.raw.crop(e),3)
            ax.set(xlim=(center[0]-20,center[0]+20),ylim=(center[1]+20,center[1]-20),aspect='equal')
            ax.add_patch(plt.Rectangle(center-13,26,26,fill=False,color='#b03090',lw=1));ax.plot(*center,'x',color='#b03090')
            ax.set_title(f'{segment} {duration}ms {r["feature"]} f{int(frame)}: {label}\n{lo-t0}..{hi-t0}us',fontsize=8)
        axes[i,1].set_xlabel(f'{r["status"]}\nscore={r["score"]:.3f}; margin={r["score_margin"]:.4f}; FB={r["forward_backward_px"]:.2f}px',fontsize=8)
    fig.suptitle('Deterministic nearest-to-seed rejection per segment/duration (all failures remain in CSV)\nA gate rejection is not evidence that a physical Marker disappeared',fontsize=12)
    fig.subplots_adjust(left=.06,right=.98,bottom=.04,top=.94,wspace=.28,hspace=.75)
    fig.savefig(OUT/'failure_cases.png',dpi=160);plt.close(fig)


def plot_actual_frames(rows):
    events,_=run.raw.load_events();times=events['t'].astype(np.int64);t0=int(times[0])
    fig,axes=plt.subplots(3,5,figsize=(20,12));colors={'F1':'#b03090','F2':'#008b45','F3':'#7b55cc'}
    for i,(segment,(start,end)) in enumerate(run.raw.SEGMENTS.items()):
        for j,frame in enumerate([0,5,10,15,19]):
            lo=t0+round(start*1e6)+frame*5000;hi=lo+5000
            e=events[(times>=lo)&(times<hi)];ax=axes[i,j];run.raw.draw_raw(ax,run.raw.crop(e),1.8)
            failed=[]
            for feature in run.SEEDS[segment]:
                r=next(r for r in rows if r['segment']==segment and r['duration_ms']==5 and r['feature']==feature and r['frame']==frame)
                if r['position_valid']:
                    center=point(r,'anchor');ax.plot(*center,'x',color=colors[feature],ms=5)
                    ax.add_patch(plt.Rectangle(center-13,26,26,fill=False,color=colors[feature],lw=.8))
                    ax.text(center[0]+1,center[1]-2,feature,color=colors[feature],fontsize=7)
                else:failed.append(feature+': lost/rejected')
            ax.set_title(f'{segment} 5ms f{frame:02d} [{(lo-t0)/1e6:.3f},{(hi-t0)/1e6:.3f})s',fontsize=9)
            ax.tick_params(labelsize=6)
            if failed:ax.set_xlabel('; '.join(failed),fontsize=7,color='#aa3333')
    fig.suptitle('Actual time-resolved review: raw events + accepted local patch anchors (x) and 27px boxes\n'
                 'No filled-in points after loss; F labels do not verify physical identity or geometric centers',fontsize=12)
    fig.subplots_adjust(left=.035,right=.995,bottom=.04,top=.90,wspace=.16,hspace=.5)
    fig.savefig(OUT/'actual_frame_review.png',dpi=180);plt.close(fig)


def report(rows,summary,dist,validation,overlap):
    lines=['# B0：局部特征观测与100ms短时关联','',
           '局部观测/短时patch关联已有实证进展：9个人工候选中7条5ms路径覆盖100ms，但2/5ms相对位移差可达9.95px。几何中心、物理身份和共同刚体仍未验证，PnP不启动。这里不是对全部约14个Marker的成功率。固定9个人工patch锚点与门槛的一次运行结果；F1/F2/F3不是MarkerID或几何中心。位置匹配只能表明局部结构在短时内满足本工程gate，不能确认同一物理Marker，更不能证明三维刚体模型。','',
           '## 数据与预先固定方法','',
           '完整[协议](protocol.md)、[人工种子与依据](seeds.json)、[布局审计](layout_audit.md)。early_turn 1.70–1.80、middle 6.30–6.40、late 11.30–11.40s。2ms各50片、5ms各20片连续不重叠；原始全序列contact sheets在运行前提供父代理审阅，冻结9种子，不替换失败样本。','',
           '参考块分别为名义中点之后[+0,+2ms)、[+0,+5ms)，块mid为+1ms/+2.5ms。两时长人工锚点同像素不表示同一瞬时物理点。27×27双极性distinct-pixel occupancy patch，sigma0.8；逐相邻块±12px整数平移搜索，cosine≥0.35、≥12不同像素、最佳与排除邻域3px后的次峰差≥0.025、前后向误差≤2px且不触搜索边。每步更新模板；任一拒绝即该方向停止，不重新接回。','',
           '累计模板更新可能漂移或换到相似结构。前后向一致性、次峰分差和非重叠预测支持都只提供自洽性，不是身份真值。原始事件质心不用于驱动匹配；几何模型不用于选点或填轨迹。匹配图像使用x450..670/y120..280的padding，以容纳搜索与patch；原始contact sheets及轨迹审阅主ROI仍为x485..630/y155..250，主ROI外支持可能未显示但在原始匹配计数中保留。','',
           '## 每个种子的完整出口','',
           '有效数包含1个人工seed；自动匹配数另列CSV。范围为有效原始块的start/end。终止原因不意味着原始结构或物理Marker消失。','',
           '|片段|F|T ms|有效/总帧|有效范围s|向前时间终止|向后时间终止|圆/椭圆有效模型数|','|---|---|---:|---:|---|---|---|---:|']
    for r in summary:
        lines.append(f'|{r["segment"]}|{r["feature"]}|{r["duration_ms"]}|{r["valid_with_manual_seed"]}/{r["frames"]}|{r["valid_start_s"]:.3f}–{r["valid_end_s"]:.3f}|{r["forward_termination"]}|{r["backward_termination"]}|{r["circle_valid"]}/{r["ellipse_valid"]}|')
    lines+=['','5ms七条局部结构可覆盖完整100ms：早F1/F2、中F1/F2、晚F1/F2/F3；早F3和中F3有歧义终止。2ms较多提前停止，不能把“2ms匹配未过gate”解释为无Marker信号。晚段5ms仍有可持续局部结构，因此上轮完整球弧失败不等于局部观测完全丢失。','',
            '拒绝统计：'+json.dumps(validation['status_counts'],ensure_ascii=False)+'。17次拒绝全部涉及歧义分差或前后向门槛，而非本轮score/最小像素/搜索边界门槛。尤其早5ms F3分差0.02469、中5ms F3分差0.02122，均接近预先固定0.025门槛；本轮保留拒绝、不事后放宽。这说明工程gate可能偏严格，不能断言其为物理遮挡。','',
            '## 2ms与5ms的位置一致性','',
            '只在5ms有效mid取相邻两片均有效的2ms轨迹线性插值，不跨失效段。原始锚点差受seed时差影响；相对位移差只减去最接近名义中点的一个共同有效时刻偏置（逐F固定，不每帧对齐）。以下为中位/P95/最大欧氏差，单位px；共同参考行的零相对差由定义产生，不能算额外验证。','',
            '|片段|F|配对数|原始锚点差 中位/P95/最大|固定初始偏置后的位移差 中位/P95/最大|','|---|---|---:|---|---|']
    for segment in run.raw.SEGMENTS:
        for feature in run.SEEDS[segment]:
            q={r['diagnostic']:r for r in dist if r['segment']==segment and r['feature']==feature and r['duration_ms']=='2_vs_5'}
            a,b=q['anchor_raw_distance'],q['anchor_relative_displacement_distance']
            lines.append(f'|{segment}|{feature}|{a["n"]}|{f(a["median"])}/{f(a["p95"])}/{f(a["max"])}|{f(b["median"])}/{f(b["p95"])}/{f(b["max"])}|')
    lines+=['','## 跨F的重叠诊断','',
            '下表只比较同时有效的人工锚点和27×27patch；IoU是正方形patch的几何交并比，不是Marker轮廓重叠率。重叠会共享事件支持，F不能直接当独立Marker约束；即使IoU为零，也不能单凭分散位置证明不同物理身份。','',
            '|片段|T|F对|同时有效帧|最小锚点距离px|最大patch IoU|','|---|---:|---|---:|---:|---:|']
    for segment in run.raw.SEGMENTS:
        for duration in [2,5]:
            for a,b in [('F1','F2'),('F1','F3'),('F2','F3')]:
                q=[r for r in overlap if r['segment']==segment and r['duration_ms']==duration and r['feature_a']==a and r['feature_b']==b]
                lines.append(f'|{segment}|{duration}|{a}/{b}|{len(q)}|{f(min(r["anchor_distance_px"] for r in q)) if q else "NA"}|{f(max(r["patch_iou"] for r in q)) if q else "NA"}|')
    lines+=['','模型与质心的跨时长比较同存cross_window_positions.csv及diagnostic_distributions.csv；没有可对齐有效中心时不填值。整数位移量化、逐帧模板更新的累积漂移、两种窗口的事件内容差异与局部ROI截断都可能贡献跨窗差；最高约10px差不能当真实Marker中心漂移或中心估计误差。轨迹差不是绝对位置误差，也不能证明共同刚体。2D投影距离可随姿态变化，本轮没有把2D距离恒定作为刚体必要条件。','',
            '## 极性、局部中心模型与时间内最大差','',
            '质心按不同像素等权计算，ON/OFF/combined分别记录；它们是不同边缘响应的空间位置统计量，不是几何中心或中心估计误差；不能通过平均ON/OFF质心就认定几何中心。圆/椭圆拟合仅对原始不同像素事件位置作候选边缘模型，至少15点、中心bounds锚点±8px、半径/半轴3..18px；至少8/12角bin、最长空角≤90°、median径向残差≤1.5px、RMS≤2.5px，椭圆轴比≤2且不触bounds。无效时正式中心为NaN，原始拟合与拒绝原因仍留在observations.csv。模型角覆盖本身依赖拟合，不能确认物理边界。','',
            '下表为有效观测里的ON−OFF质心距离和圆−椭圆中心距离，中位/P95/最大；后一列只比较双方模型gate有效帧。5ms各候选ON/OFF质心距离中位为5.43–16.06px，反映不同极性边缘响应。圆/椭圆同源拟合差小并非独立准确性证据，亦不能排除共同偏差或结构混合。早F2的5ms全部20块可匹配但圆0/20、椭圆1/20有效，尤其说明可匹配patch不等于可用几何中心。','',
            '|片段|F|T|ON−OFF n / 中位/P95/最大|圆−椭圆 n / 中位/P95/最大|','|---|---|---:|---|---|']
    for segment in run.raw.SEGMENTS:
        for feature in run.SEEDS[segment]:
            for duration in [2,5]:
                q={r['diagnostic']:r for r in dist if r['segment']==segment and r['feature']==feature and r['duration_ms']==duration}
                a=q.get('on_centroid_minus_off_centroid',dict(n=0,median=np.nan,p95=np.nan,max=np.nan));b=q.get('circle_minus_ellipse',dict(n=0,median=np.nan,p95=np.nan,max=np.nan))
                lines.append(f'|{segment}|{feature}|{duration}|{a["n"]} / {f(a["median"])}/{f(a["p95"])}/{f(a["max"])}|{b["n"]} / {f(b["median"])}/{f(b["p95"])}/{f(b["max"])}|')
    lines+=['','diagnostic_distributions.csv还报告各模型/质心相对锚点或彼此的位置差，以及同轨迹这些偏移在时间内的最大两两变化直径offset_change_diameter_px；不是物理中心漂移真值。观测表future_*字段实际指传播方向下一片非重叠块：direction=+1为时间未来，direction=−1为更早；协议中“未来”应按此实现术语澄清。它只在常速预测patch上评分，不搜索、不用于当前接受；加速、形变和遮挡均可降低它，错误持续结构也可升高。','',
            '## 阶段出口','',
            '- 局部结构可持续追踪：上述七条5ms完整路径仅在人工种子、±12px搜索和当前门槛下成立；两条5ms压力候选及多数2ms路径提前拒绝。不可把缺失补成完整Marker轨迹。',
            '- 几何中心未确认：patch锚点、极性质心、局部圆/椭圆有效模型是不同观测定义；即使模型通过角覆盖gate，也不证明它是球形Marker的投影中心。单独可见弧、重叠事件和未来形变会改变位置。',
            '- 物理身份未验证：F名仅是局部结构标签，累计更新可能跟随边缘而非同一球面点。通过FB与跨窗对齐不证明物理身份。',
            '- 不具备宣称姿态的条件：见[布局审计](layout_audit.md)，已找到10点r20mm模板但未绑定当前约14点球，当前G1相机K亦未找到；本轮不做PnP、三维刚体或角速度结论。','',
            '## 复现和独立检查','',
            f'630观察行唯一且完整；{validation["independently_rescored_pairs"]}个代表接受/拒绝匹配用绝对时间boolean选窗、显式像素集合累积及逐位移直接cosine独立核对最佳与次峰，均匹配。36个方向逐条验证丢失后未重连。人工已知(4,3)px位移恢复、相似双峰歧义拒绝、空patch拒绝、重现后不接回均通过。原始210个窗的独立时间检查在preview_metadata.json。完整原始来源/配置/协议/种子/代码hash在metadata.json和validation.json。','',
            '```sh','g1_spin/.venv/bin/python g1_spin/experiments/marker_b0_20261101/continuous_preview.py','g1_spin/.venv/bin/python g1_spin/experiments/marker_b0_20261101/run.py','g1_spin/.venv/bin/python g1_spin/experiments/marker_b0_20261101/summarize.py','```','',
            '图：continuous_*_{2,5}ms.png为全部连续原始片；trajectories.png为人工锚点路径；actual_frame_review.png为每段5ms第0/5/10/15/19帧的实际事件和当帧有效patch（失效不补点）；cross_window_differences.png为对齐诊断；failure_cases.png按每片段/时长距seed最近拒绝确定性展示，全部17次拒绝和250帧停止后状态在observations.csv。无真实中心/身份标签，无绝对误差、精确率/召回率、三维旋转结论。本轮仅写此新目录，未改旧产物，未安装依赖、训练、联网或操作硬件。']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':main()
