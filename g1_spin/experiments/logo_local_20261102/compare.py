"""B2 joint availability only: independent tracks remain unchanged, no multi-point estimator."""
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import review
import exploratory
from validate_review import read_csv

OUT=review.OUT
MANUAL={
 ('early_turn','F1'):('round_candidate','双极性圆样边缘可见；27px框截断/内部填充使圆模型支持随时间变化'),
 ('early_turn','F2'):('ambiguous','混合多弧、初始主要位于框外；不可由不圆或拟合失败推商标'),
 ('early_turn','F3'):('round_partial_candidate','起始圆/椭圆边缘，随后变窄和离开ROI；丢失后不作同一结构解释'),
 ('middle','F1'):('round_candidate','原始圆样边缘连续，部分窗口圆/椭圆模型支持'),
 ('middle','F2'):('round_candidate','圆样双极性结构，锚点框常截去一侧；模型少不代表logo'),
 ('middle','F3'):('ambiguous','薄弧/稀疏极性结构，不能排除斜视圆点或单一边缘'),
 ('late','F1'):('round_candidate','分散的双极性圆/椭圆边缘，多数局部模型gate有效'),
 ('late','F2'):('round_candidate','分散圆/椭圆边缘，多数模型gate有效'),
 ('late','F3'):('round_partial_candidate','近下方截断的圆/窄弧候选，边界/遮挡风险保留'),
 ('post_survey_11_75','L'):('irregular_logo_candidate_UNCONFIRMED','全片survey后选取；多个连续原始块有非单薄弧的多方向细节，符合用户形态提示；邻近圆弧/重叠尚未排除'),
 ('post_survey_11_75','A1'):('round_candidate','与L约53px的较远圆弧对照，不称附近'),
 ('post_survey_11_75','A2'):('round_partial_candidate','与L约64px的较远窄弧/下方截断对照'),
 ('post_survey_11_75','N1'):('round_candidate','与L约23px的附近圆/椭圆候选，27px支持可混入L'),
 ('post_survey_11_75','N2'):('round_partial_candidate','与L约26px的附近窄弧候选，可能共享局部事件')}


def usable(r):return r['position_valid'] is True

def geometric(r):return usable(r) and (r.get('circle_valid') is True or r.get('ellipse_valid') is True)

def number_stats(values):
    values=np.array(values,dtype=float);values=values[np.isfinite(values)]
    return (float(np.median(values)),float(np.quantile(values,.95)),float(np.max(values))) if len(values) else (np.nan,np.nan,np.nan)

def f(x):return 'NA' if not np.isfinite(x) else f'{x:.3f}'


def main():
    old=review.read_rows();new=read_csv('exploratory_observations.csv')+read_csv('near_observations.csv')
    evidence=[];tracks=[]
    for segment,feature in MANUAL:
        source=old if segment!='post_survey_11_75' else new
        q=[r for r in source if r['segment']==segment and r['feature']==feature and r['duration_ms']==5]
        valid=[r for r in q if usable(r)];cls,reason=MANUAL[segment,feature]
        evidence.append(dict(segment=segment,feature=feature,human_candidate_level=cls,human_basis=reason,physical_identity_verified=False,
                             selection='old manually conditioned B0 sample' if source is old else 'post-survey exploratory; NOT independent test',
                             morphology_scope_ms=100 if source is old else 50,total_5ms_frames=len(q),valid_5ms=len(valid),
                             circle_gate_valid=sum(r.get('circle_valid') is True for r in valid),ellipse_gate_valid=sum(r.get('ellipse_valid') is True for r in valid),
                             circle_or_ellipse_gate_valid=sum(geometric(r) for r in valid),
                             median_border_2px_fraction=float(np.median([r['boundary_2px_pixel_fraction'] for r in valid])) if source is not old else old_median(segment,feature,'boundary_2px_pixel_fraction'),
                             descriptor_is_not_classifier=True))
        for duration in [5,2]:
            qq=[r for r in source if r['segment']==segment and r['feature']==feature and r['duration_ms']==duration and r['frame']>=(10 if duration==5 else 25)]
            vv=[r for r in qq if usable(r)];rejections=[r for r in qq if str(r['status']).startswith('rejected_')]
            tracks.append(dict(segment=segment,feature=feature,duration_ms=duration,evaluation_ms=50,frames=len(qq),valid_including_seed=len(vv),
                               circle_valid=sum(r.get('circle_valid') is True for r in vv),ellipse_valid=sum(r.get('ellipse_valid') is True for r in vv),
                               either_model_valid=sum(geometric(r) for r in vv),first_rejected_frame=rejections[0]['frame'] if rejections else np.nan,
                               rejection_reason=rejections[0]['status'] if rejections else 'none',
                               minimum_accepted_score_margin=min((r.get('score_margin',np.inf) for r in vv if r['status']=='accepted'),default=np.nan)))
    review.write_csv('candidate_evidence.csv',evidence);review.write_csv('forward_track_summary.csv',tracks)
    definitions=[('old_ambiguous_early','early_turn','F2',['F1','F3'],old),('old_ambiguous_middle','middle','F3',['F1','F2'],old),
                 ('new_distant_control','post_survey_11_75','L',['A1','A2'],new),('new_NEAR','post_survey_11_75','L',['N1','N2'],new)]
    joint=[];summary=[]
    for branch,segment,primary,auxiliaries,source in definitions:
        for duration in [5,2]:
            lookup={(r['feature'],int(r['frame'])):r for r in source if r['segment']==segment and r['duration_ms']==duration}
            frames=range(10,20) if duration==5 else range(25,50)
            for n in [0,1,2]:
                features=[primary]+auxiliaries[:n];rr=[]
                for k in frames:
                    p=lookup[primary,k];others=[lookup[a,k] for a in auxiliaries[:n]]
                    passed=usable(p) and all(usable(r) for r in others)
                    model_pass=usable(p) and all(geometric(r) for r in others) if others else np.nan
                    record=dict(branch=branch,segment=segment,duration_ms=duration,frame=k,midpoint_s=p['midpoint_s'],group='+'.join(features),auxiliary_count=n,
                                primary_human_level=MANUAL[segment,primary][0],primary_status=p['status'],primary_position_valid=usable(p),
                                primary_anchor_x=p['anchor_x'],primary_anchor_y=p['anchor_y'],joint_position_available=passed,joint_auxiliary_models_available=model_pass,
                                missing_position_features=','.join(r['feature'] for r in [p]+others if not usable(r)),
                                missing_auxiliary_model_features=','.join(r['feature'] for r in others if not geometric(r)),
                                primary_output_changed_by_auxiliaries=False)
                    joint.append(record);rr.append(record)
                summary.append(dict(branch=branch,duration_ms=duration,group='+'.join(features),total_frames=len(rr),
                                    primary_available=sum(r['primary_position_valid'] for r in rr),joint_position_available=sum(r['joint_position_available'] for r in rr),
                                    joint_auxiliary_models_available=sum(r['joint_auxiliary_models_available'] is True for r in rr) if n else np.nan,
                                    coverage_loss_while_primary_valid=sum(r['primary_position_valid'] and not r['joint_position_available'] for r in rr),
                                    primary_matching_changed=False,scope='50ms forward; joint condition only, not a joint estimator'))
    review.write_csv('joint_observations.csv',joint);review.write_csv('joint_summary.csv',summary)
    frames,checks=exploratory.windows();overlap=[]
    for duration in [5,2]:
        lookup={(r['feature'],int(r['frame'])):r for r in new if r['duration_ms']==duration}
        for k in (range(10,20) if duration==5 else range(25,50)):
            p=lookup['L',k]
            if not usable(p):continue
            for aux in ['A1','A2','N1','N2']:
                a=lookup[aux,k]
                if not usable(a):continue
                pc=np.array([p['anchor_x'],p['anchor_y']]);ac=np.array([a['anchor_x'],a['anchor_y']]);delta=abs(pc-ac)
                inter=float(np.prod(np.maximum(0,27-delta)));e=frames[duration][k]['events']
                pm=(abs(e['x'].astype(float)-pc[0])<=13)&(abs(e['y'].astype(float)-pc[1])<=13)
                am=(abs(e['x'].astype(float)-ac[0])<=13)&(abs(e['y'].astype(float)-ac[1])<=13)
                xy=lambda z:set(zip(z['x'].tolist(),z['y'].tolist()))
                shared=len(xy(e[pm])&xy(e[am]))
                overlap.append(dict(duration_ms=duration,frame=k,auxiliary=aux,branch='near' if aux.startswith('N') else 'distant',
                                    anchor_distance_px=float(np.linalg.norm(delta)),square27_patch_iou=inter/(1458-inter),
                                    shared_raw_events=int(np.sum(pm&am)),shared_distinct_xy_pixels=shared,
                                    shared_fraction_of_L_events=float(np.sum(pm&am)/np.sum(pm)) if np.sum(pm) else np.nan))
    review.write_csv('patch_overlap.csv',overlap)
    aligned=align_feature(new,'L')+align_feature(new,'N1');review.write_csv('cross_duration_positions.csv',aligned)
    review.write_csv('L_cross_duration.csv',[r for r in aligned if r['feature']=='L'])
    validation=validate(old,new,joint,frames)
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    metadata=dict(user_rule='Irregular appearance is logo, round is attached Marker: treated as human candidate hint only',
                  no_identity_truth=True,no_automatic_semantic_classifier=True,original_nine_not_replaced=True,
                  primary_source_sha256=hashlib.sha256((review.B0/'observations.csv').read_bytes()).hexdigest(),
                  review_metadata_sha256=hashlib.sha256((OUT/'review_metadata.json').read_bytes()).hexdigest(),
                  exploratory_protocol_sha256=hashlib.sha256((OUT/'exploratory_protocol.json').read_bytes()).hexdigest(),
                  near_protocol_sha256=hashlib.sha256((OUT/'near_protocol.json').read_bytes()).hexdigest(),
                  compare_sha256=hashlib.sha256((OUT/'compare.py').read_bytes()).hexdigest(),
                  survey_note='Parent fixed 35-window survey selected new 11.75s candidate; post-selection exploratory addition, no independent test')
    (OUT/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    plot_joint(joint);plot_cross_duration(aligned);report(evidence,tracks,summary,overlap,aligned,validation)
    print(json.dumps({k:v for k,v in validation.items() if k!='score_checks'},indent=2))


def old_median(segment,feature,field):
    q=[r[field] for r in read_csv('morphology.csv') if r['segment']==segment and r['feature']==feature and r['duration_ms']==5 and usable(r)]
    return float(np.median(q))


def align_feature(rows,feature):
    a=[r for r in rows if r['feature']==feature and r['duration_ms']==2];b=[r for r in rows if r['feature']==feature and r['duration_ms']==5]
    times=np.array([r['midpoint_s'] for r in a]);result=[]
    for r in b:
        if not usable(r):continue
        i=int(np.searchsorted(times,r['midpoint_s']))
        if not 0<i<len(a) or not usable(a[i-1]) or not usable(a[i]):continue
        w=(r['midpoint_s']-times[i-1])/(times[i]-times[i-1]);p=(1-w)*np.array([a[i-1]['anchor_x'],a[i-1]['anchor_y']])+w*np.array([a[i]['anchor_x'],a[i]['anchor_y']])
        delta=p-[r['anchor_x'],r['anchor_y']]
        result.append(dict(feature=feature,midpoint_s=r['midpoint_s'],dx=delta[0],dy=delta[1],raw_distance=float(np.linalg.norm(delta))))
    if result:
        bias=np.array([result[0]['dx'],result[0]['dy']])
        for r in result:r.update(common_reference_s=result[0]['midpoint_s'],relative_displacement_distance=float(np.linalg.norm(np.array([r['dx'],r['dy']])-bias)))
    return result


def validate(old,new,joint,frames):
    meta=json.loads((OUT/'review_metadata.json').read_text());near=json.loads((OUT/'near_metadata.json').read_text())
    assert review.digest(review.B0)==meta['previous_B0_sha256']
    assert hashlib.sha256((OUT/'exploratory_observations.csv').read_bytes()).hexdigest()==near['existing_exploratory_observations_sha256']
    assert len(new)==175==len({(r['duration_ms'],r['feature'],r['frame']) for r in new})
    assert all(not r['primary_output_changed_by_auxiliaries'] for r in joint)
    events,_=review.b0.raw.load_events();times=events['t'].astype(np.int64);t0=int(times[0]);generalized_checks=0
    for segment,(start,end) in review.b0.raw.SEGMENTS.items():
        for duration in [5,2]:
            original=[]
            for k in range(100//duration):
                lo=t0+round(start*1e6)+k*duration*1000;hi=lo+duration*1000
                original.append(dict(events=events[(times>=lo)&(times<hi)],start_s=(lo-t0)/1e6,end_s=(hi-t0)/1e6))
            images=[review.b0.image(f['events']) for f in original]
            for feature,seed in review.b0.SEEDS[segment].items():
                rr=exploratory.forward_general(original,images,seed,duration,feature)
                prior={r['frame']:r for r in old if r['segment']==segment and r['duration_ms']==duration and r['feature']==feature}
                for r in rr:
                    p=prior[r['frame']];assert r['status']==p['status']
                    for key in ['anchor_x','anchor_y','score','second_score','score_margin','forward_backward_px']:
                        if key in r and np.isfinite(r[key]):assert r[key]==p[key]
                    generalized_checks+=1
    score_checks=[]
    for duration in [5,2]:
        for feature in ['L','A1','A2','N1','N2']:
            q=[r for r in new if r['duration_ms']==duration and r['feature']==feature];lost=False
            for r in q:
                if lost:assert r['status']=='not_attempted_after_loss' and not usable(r)
                if str(r['status']).startswith('rejected_'):lost=True
            selected=[]
            for pattern in ['accepted','rejected_']:
                candidate=next((r for r in q if str(r['status']).startswith(pattern)),None)
                if candidate:selected.append(candidate)
            # Independent occupancy accumulation via tuple sets, then finite-radius smoothing.
            def image(e):
                im=np.zeros((2,161,221))
                for x,y,p in set(zip(e['x'].tolist(),e['y'].tolist(),(e['p']>0).tolist())):
                    if 450<=x<=670 and 120<=y<=280:im[int(p),y-120,x-450]=1
                return np.stack([__import__('scipy').ndimage.gaussian_filter(channel,.8,radius=3) for channel in im])
            for r in selected:
                previous=next(x for x in q if x['frame']==r['source_frame']);center=np.array([previous['anchor_x'],previous['anchor_y']])
                source=image(frames[duration][int(r['source_frame'])]['events']);target=image(frames[duration][int(r['frame'])]['events'])
                template=review.b0.patch(source,center);scores=np.zeros((25,25))
                for iy,dy in enumerate(range(-12,13)):
                    for ix,dx in enumerate(range(-12,13)):
                        p=review.b0.patch(target,center+[dx,dy]);den=np.linalg.norm(template)*np.linalg.norm(p);scores[iy,ix]=np.sum(template*p)/den if den else 0
                iy,ix=np.unravel_index(np.argmax(scores),scores.shape)
                assert np.array_equal(center+[ix-12,iy-12],[r['proposed_x'],r['proposed_y']]) and abs(scores[iy,ix]-r['score'])<1e-10
                scores[max(0,iy-3):iy+4,max(0,ix-3):ix+4]=-np.inf
                assert abs(np.max(scores)-r['second_score'])<1e-10
                score_checks.append(dict(duration_ms=duration,feature=feature,frame=r['frame'],score_map_and_second_match=True))
    return dict(original_B0_files_unchanged=True,original_L_A1_A2_observations_unchanged=True,
                generalized_forward_loop_matches_original_B0_rows=generalized_checks,new_observation_rows_complete=175,
                independent_direct_score_maps=len(score_checks),score_checks=score_checks,
                no_reconnection_all_10_paths=True,joint_conditions_never_modify_primary=True,
                candidate_semantics_validated=False,post_survey_selection_not_independent_test=True)


def plot_joint(rows):
    fig,axes=plt.subplots(2,2,figsize=(13,7))
    for i,duration in enumerate([5,2]):
        for j,branch in enumerate(['new_distant_control','new_NEAR']):
            ax=axes[i,j]
            for n in [0,1,2]:
                q=[r for r in rows if r['branch']==branch and r['duration_ms']==duration and r['auxiliary_count']==n]
                ax.plot([r['midpoint_s'] for r in q],[n+(0.7 if r['joint_position_available'] else 0) for r in q],'o-',ms=3,label=q[0]['group'])
            ax.set(yticks=[.35,1.35,2.35],yticklabels=['single L','L + 1 aux','L + 2 aux'],title=f'{branch} {duration}ms (raised = joint positions usable)')
            ax.grid(alpha=.2)
    fig.suptitle('Joint quality conditions change coverage; they do NOT change or refill L matching',fontsize=12)
    fig.tight_layout(rect=(0,0,1,.94));fig.savefig(OUT/'joint_coverage.png',dpi=180);plt.close(fig)


def plot_cross_duration(rows):
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for ax,feature in zip(axes,['L','N1']):
        q=[r for r in rows if r['feature']==feature]
        ax.plot([r['midpoint_s'] for r in q],[r['raw_distance'] for r in q],'o-',label='raw anchor difference')
        ax.plot([r['midpoint_s'] for r in q],[r['relative_displacement_distance'] for r in q],'s-',label='after ONE fixed initial offset')
        ax.set(title=feature+' 2ms vs 5ms, valid adjacent interpolation only',xlabel='Aligned time (s)',ylabel='Distance (px)');ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle('Full tracking availability does not establish stable geometric position; no identity/center truth',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.91));fig.savefig(OUT/'cross_duration_differences.png',dpi=180);plt.close(fig)


def report(evidence,tracks,summary,overlap,aligned,validation):
    lines=['# B2：不规则商标候选与少量圆样辅助观测','',
           '旧9个人工种子仍不能可靠指认商标；将circle失败或薄弧叫logo会误判。全片固定survey后，在11.75s发现具有连续非单薄弧细节的不规则商标候选L。新片段为post-survey探索性增补，符合用户形态提示，尚未排除圆点斜视、重叠和ROI截断，不能称物理商标身份已确认。','',
           'L在11.750–11.800s的5ms共10/10块、2ms共25/25块通过B0原匹配gate（均含1人工seed），是局部形态/关联延续证据，不是商标识别准确率。然而减一次初始偏置后L的2/5ms相对位移差仍最大9.520px，N1最大7.977px，尚未达到稳定几何位置观测。辅助点仍独立跟踪，既不改变也不回填L；联合门槛只改变整组观测覆盖，没有证明中心、姿态或轴准确性提高。','',
           '## 人工提示、原始证据与选择过程','',
           '用户提示“不规则外观=商标、圆形=贴点”用于人工候选分层；并非以任意形态阈值训练分类器。旧主样本固定不替换。先审[参考全ROI](reference_context.png)和card_*局部100ms演化卡：上排5ms帧0/5/10/15/19，下排最近2ms原始窗，mid差±0.5ms。失效时明确固定seed审阅，绝不补成同一特征轨迹。','',
           '父代理按全片0.25+0.5k的35个固定5ms窗[定时survey](whole_recording_survey_1.png)（另whole_recording_survey_2.png）检索，随后审[三个连续窗细节](irregular_survey_details.png)，再提出11.70–11.80s探索分支。新种子L(582,212)、较远A1(533,192)、A2(523,238)先于新跟踪冻结；之后父代理仅看原始seed图、尚未见跟踪结果时提出near N1(588,190)、N2(576,237)，单独追加near_protocol.json，原三条输出逐字保持。新数据选择不是独立测试，不能据此估计泛化率。','',
           '先看[完整5ms原始序列](exploratory_raw_5ms.png)、[完整2ms对照](exploratory_raw_2ms.png)、[跟踪前原始候选卡](exploratory_seed_review.png)，再看[较远对照跟踪](exploratory_tracking_review.png)与[near跟踪](near_tracking_review.png)。L框有多方向细节，并非仅单一细弧，但上方圆弧/右侧相邻结构仍可能混入。','',
           '## 候选证据卡：人工推断与几何gate分开','',
           '旧样本下表为100ms的20个5ms块；新分支为向前50ms的10块。circle/ellipse数量仅是原B0局部模型可拟合gate，不能验证语义；连通分量是事件碎片，不是Marker数；边框占比描述27px截断风险，不是物理遮挡真值。','',
           '|片段/F|人工候选层级|有效块/总块|circle/ellipse gate有效|人工依据/替代解释|','|---|---|---:|---:|---|']
    for r in evidence:lines.append(f'|{r["segment"]}/{r["feature"]}|{r["human_candidate_level"]}|{r["valid_5ms"]}/{r["total_5ms_frames"]}|{r["circle_gate_valid"]}/{r["ellipse_gate_valid"]}|{r["human_basis"]}|')
    lines+=['','所有原始形态指标及失败保留在morphology.csv；候选层级与人工依据在candidate_evidence.csv。早F2混合/截断、middle F3薄弧均保持ambiguous；middle F2圆样边缘常被27px框截掉一侧，说明模型失败不能推出logo。新L的圆/椭圆gate均失败也不是判定商标的依据。','',
            '## 新分支固定门槛前向50ms结果','',
            '仍为27×27双极性不同像素occupancy、sigma0.8、±12px整数平移、≥12像素、cosine≥0.35、次峰（排除3px邻域）分差≥0.025、FB≤2px、不触搜索边。参考块5ms为[11.750,11.755)，2ms为[11.750,11.752)，不是同一瞬间；一次拒绝后不接回。没有固定贴点布局、共同平移或平面假设。','',
            '|特征|T ms|有效/总块|circle/ellipse gate有效|首次拒绝|','|---|---:|---:|---:|---|']
    for r in tracks:
        if r['segment']=='post_survey_11_75':lines.append(f'|{r["feature"]}|{r["duration_ms"]}|{r["valid_including_seed"]}/{r["frames"]}|{r["circle_valid"]}/{r["ellipse_valid"]}|{r["rejection_reason"]}|')
    lines+=['','L和N1的2ms位置只在相邻有效块间线性插值到5ms mid；相对位移差仅减最接近名义参考时刻的一个固定初始偏置，绝不逐帧重对齐、不跨失效插值。下面为中位/P95/最大欧氏差，单位px；初始零相对差由定义产生，不是额外真值。','',
            '|特征|配对数|原始锚点差 中位/P95/最大|减单个初始偏置后的位移差 中位/P95/最大|','|---|---:|---|---|']
    for feature in ['L','N1']:
        q=[r for r in aligned if r['feature']==feature];raw=number_stats([r['raw_distance'] for r in q]);rel=number_stats([r['relative_displacement_distance'] for r in q])
        lines.append(f'|{feature}|{len(q)}|{f(raw[0])}/{f(raw[1])}/{f(raw[2])}|{f(rel[0])}/{f(rel[1])}/{f(rel[2])}|')
    lines+=['','整数位移量化、模板更新漂移、窗口内容和ROI截断均可贡献差异，不能解释为商标/贴点中心误差；完整位置/初始参考时间见cross_duration_positions.csv，变化图见[cross_duration_differences.png](cross_duration_differences.png)。','',
            '## 公平比较：单候选与联合可用性','',
            '每组使用完全相同主体轨迹，不平均不同点的2D位置/速度，不要求投影距离恒定，不用未验证辅助补主体。position联合gate要求指定所有轨迹当帧可用；aux-model联合gate另外要求每个辅助至少一种circle/ellipse gate有效，主体L不强套圆模型。单主体的辅助模型列为NA，不代表L几何中心已有效。','',
            '|分支|T|组|主体可用|联合position可用/总块|联合aux-model可用|主体有效时因辅助失去的覆盖块|','|---|---:|---|---:|---:|---:|---:|']
    for r in [x for x in summary if x['branch'].startswith('new')]:lines.append(f'|{r["branch"]}|{r["duration_ms"]}|{r["group"]}|{r["primary_available"]}|{r["joint_position_available"]}/{r["total_frames"]}|{int(r["joint_auxiliary_models_available"]) if np.isfinite(r["joint_auxiliary_models_available"]) else "NA"}|{r["coverage_loss_while_primary_valid"]}|')
    lines+=['','旧主样本条件组压缩如下（5ms向前50ms；2ms逐块结果留CSV），主体均为ambiguous假定不规则，不能称商标真值。','',
            '|旧条件组|单主体position|+1辅助position / aux-model|+2辅助position / aux-model|','|---|---:|---:|---:|']
    for branch in ['old_ambiguous_early','old_ambiguous_middle']:
        q=[r for r in summary if r['branch']==branch and r['duration_ms']==5]
        values=[]
        for r in q:values.append(str(r['joint_position_available']) if r['group'].count('+')==0 else f"{r['joint_position_available']} / {int(r['joint_auxiliary_models_available'])}")
        lines.append(f'|{branch}|{values[0]}/10|{values[1]}|{values[2]}|')
    lines+=['','旧early/middle两组主体仍是ambiguous假定不规则候选，不能称商标对照真值。新near组更贴近“商标附近点”的构想；A1/A2距L约53/64px，仅为较远分散round对照，N1/N2约23/26px。辅助轨迹提供其自身事件/形态延续和局部模型质量证据，不自动降低L匹配歧义。联合gate是逐帧质量标记，可因模型gate波动重新有效；匹配真正丢失的ID始终不复活。','',
            '5ms near N1/N2虽均10/10块可跟踪，但其局部圆/椭圆模型均0/10有效；同时有两条辅助轨迹不等于增加两个可信几何中心。patch模板更新可能沿边缘游移，尤其N2锚点不等于可见圆/弧的中心。2ms辅助模型有部分支持，但N2较早因歧义丢失，联合覆盖相应下降。','',
            '## near支持重叠：不能称独立几何约束','',
            '|辅助|T|同时有效块|最小锚点距离px|最大27px patch IoU|最多共享原始事件/不同像素|最多共享L事件比例|','|---|---:|---:|---:|---:|---:|---:|']
    for aux in ['N1','N2','A1','A2']:
        for duration in [5,2]:
            q=[r for r in overlap if r['auxiliary']==aux and r['duration_ms']==duration]
            lines.append(f'|{aux}|{duration}|{len(q)}|{f(min(r["anchor_distance_px"] for r in q))}|{f(max(r["square27_patch_iou"] for r in q))}|{max(r["shared_raw_events"] for r in q)}/{max(r["shared_distinct_xy_pixels"] for r in q)}|{max(r["shared_fraction_of_L_events"] for r in q):.1%}|')
    lines+=['','共享支持直接由同一原始窗的事件索引/像素交集计数，说明邻近patch不能直接当独立轴约束。IoU为零也不证明不同物理身份或独立统计噪声。本轮没有由多轨迹数量推导三维约束或识别准确性。','',
            '## 仍不能确定与下一步','',
            '用户形态提示在新L原始事件上有可审阅落实：稳定不规则候选可按原门槛短时跟踪；旧9样本不能可靠确认商标这一负结果保留。尚未确认L确为商标而非混合/投影结构，也未确认patch锚点是同一物理球面点或几何中心。','',
            '目前不能宣称辅助改善主体匹配：方法仍是独立跟踪，联合筛选会损失覆盖，附近patch共享事件且5ms局部圆中心支持不足。辅助轨迹没有用于修正L，因此本轮也不能判断经验证观测模型下的辅助融合能否提升精度。需要人工跨帧确认商标/圆点边界与重叠、真实遮挡/出入视野及不同布局记录，才可评价身份与失效检测；下一步先分离整数位移量化、递推模板累积漂移与ROI混入的贡献：固定参考外观模板对递推模板、亚像素对整数位移、不同ROI支持各作独立对照，不同时优化。再检查L内部可重复角点/笔画交汇及辅助Marker独立边界，验证何者对应同一物理位置；不急继续加点，不能将整个patch锚点直接提升到球面或送PnP。需要经验证的观测模型才能利用辅助构造刚体或姿态约束。PnP、三维旋转与轴准确性均未启动。','',
            '## 复现与核验','',
            f'旧630行、90个审阅panel、210个绝对时间窗保留；363个有效形态行由原始xy集合和独立8邻域flood-fill/边框坐标核对。新175行全保留；适用于2/5ms长度的前向循环实现与原B0的{validation["generalized_forward_loop_matches_original_B0_rows"]}行完全相同；{validation["independent_direct_score_maps"]}个新接受/拒绝匹配25×25分数图与次峰直接内积复算通过。L/A1/A2原三条结果在near追加后hash未变，全部10路径丢失后未重连，联合条件不改主体。语义身份没有被这些实现检查验证。','',
            '```sh','g1_spin/.venv/bin/python g1_spin/experiments/logo_local_20261102/review.py','g1_spin/.venv/bin/python g1_spin/experiments/logo_local_20261102/validate_review.py','g1_spin/.venv/bin/python g1_spin/experiments/logo_local_20261102/exploratory.py','g1_spin/.venv/bin/python g1_spin/experiments/logo_local_20261102/exploratory.py --track','g1_spin/.venv/bin/python g1_spin/experiments/logo_local_20261102/nearby.py','g1_spin/.venv/bin/python g1_spin/experiments/logo_local_20261102/compare.py','```','',
            '父代理survey复现脚本whole_recording_survey.py与survey_followup.py另存本目录。协议为protocol.md、exploratory_protocol.json、near_protocol.json；原始来源part hash见review_metadata，脚本/配置hash见metadata与validation及各分支metadata。未改原始数据/B0/B1/旋转，未训练、联网、装依赖或操作硬件。']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':main()
