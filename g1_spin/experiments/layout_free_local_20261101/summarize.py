"""Independent payload, lookup-context and pressure-outcome checks for fixed B1 inputs."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter
import run

OUT=Path(__file__).resolve().parent


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


def f(x):return 'NA' if not np.isfinite(x) else f'{x:.3f}'


def main():
    baseline=read('baseline.csv');rows=read('observations.csv');interventions=read('interventions.csv');contexts=read('context_checks.csv')
    summaries=[]
    for segment in run.b0.raw.SEGMENTS:
        for variant in run.VARIANTS:
            q=[r for r in rows if r['segment']==segment and r['variant']==variant]
            accepted=[r for r in q if r['position_valid']];rejected=[r for r in q if str(r['status']).startswith('rejected_')]
            payload=[r for r in interventions if r['segment']==segment and r['variant']==variant and r['intervention_active']]
            copy_dist=[r['accepted_distance_from_copy_anchor'] for r in accepted if np.isfinite(r.get('accepted_distance_from_copy_anchor',np.nan))]
            summaries.append(dict(segment=segment,variant=variant,evaluation_ms=50,frames=10,valid_including_manual_seed=len(accepted),
                                  first_rejected_frame=rejected[0]['frame'] if rejected else np.nan,
                                  rejection_reason=rejected[0]['status'] if rejected else 'none',
                                  not_attempted_after_loss=sum(r['status']=='not_attempted_after_loss' for r in q),
                                  max_accepted_deviation_from_baseline=max(r.get('accepted_distance_from_baseline',0) for r in accepted),
                                  min_accepted_distance_to_copy_anchor=min(copy_dist) if copy_dist else np.nan,
                                  copy_anchor_attraction_count=sum(r.get('copy_anchor_attraction',False) is True for r in q),
                                  total_added_events=sum(r['added_events'] for r in payload),total_removed_events=sum(r['removed_events'] for r in payload),
                                  changed_occupancy_min=min(r['occupancy_changed_polarity_pixels'] for r in payload),
                                  changed_occupancy_max=max(r['occupancy_changed_polarity_pixels'] for r in payload)))
    run.write_csv('summary.csv',summaries)
    validation,neighbor=validate(baseline,rows,interventions,contexts)
    run.write_csv('neighbor_patch_counts.csv',neighbor)
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    report(baseline,rows,interventions,summaries,validation,neighbor)
    print(json.dumps({k:v for k,v in validation.items() if k not in ['independent_scores','payload_checks']},indent=2))


def validate(baseline,rows,interventions,contexts):
    meta=json.loads((OUT/'metadata.json').read_text());events,manifest=run.b0.raw.load_events();times=events['t'].astype(np.int64);t0=int(times[0]);old=run.read_b0()
    assert len(baseline)==90==len({(r['segment'],r['feature'],r['frame']) for r in baseline})
    assert len(rows)==90==len({(r['segment'],r['variant'],r['frame']) for r in rows})
    assert run.digest(run.B0)==meta['baseline_source_files_sha256']
    assert hashlib.sha256((OUT/'run.py').read_bytes()).hexdigest()==meta['script_sha256']
    assert hashlib.sha256((OUT/'protocol.md').read_bytes()).hexdigest()==meta['protocol_sha256']
    for r in baseline:
        prior=old[r['segment'],r['feature'],int(r['frame'])]
        assert r['status']==prior['status']
        for key in ['anchor_x','anchor_y','proposed_x','proposed_y','score','second_score','score_margin','forward_backward_px','reverse_score']:
            if key in r and np.isfinite(r[key]):assert r[key]==float(prior[key])
    payload_checks=[];neighbor=[];independent_scores=[]
    for segment,(start,end) in run.b0.raw.SEGMENTS.items():
        originals={}
        for k in range(10,20):
            lo=t0+round(start*1e6)+k*5000;hi=lo+5000
            originals[k]=events[(times>=lo)&(times<hi)]
        for variant in run.VARIANTS:
            path=OUT/f'{segment}_{variant}_events.npz'
            assert hashlib.sha256(path.read_bytes()).hexdigest()==meta['intervention_archives_sha256'][path.name]
            with np.load(path) as data:
                for k in range(10,20):
                    original=originals[k];modified=data[f'frame_{k:02d}'];added=data[f'added_{k:02d}'];removed=data[f'removed_{k:02d}']
                    center=np.array([old[segment,'F1',k]['anchor_x'],old[segment,'F1',k]['anchor_y']])
                    # Independent direct boolean definitions, not the intervention helper.
                    near=(original['x'].astype(float)>=center[0]-13)&(original['x'].astype(float)<=center[0]+13)&(original['y'].astype(float)>=center[1]-13)&(original['y'].astype(float)<=center[1]+13)
                    if variant=='occlusion' and k in run.ACTIVE:
                        assert np.array_equal(removed,original[near]) and np.array_equal(modified,original[~near])
                    elif variant=='clone_plus8' and k in run.ACTIVE:
                        source=data[f'clone_source_{k:02d}'];assert np.array_equal(source,original[near])
                        for field in original.dtype.names:
                            if field=='x':assert np.array_equal(added[field].astype(int),source[field].astype(int)+8)
                            else:assert np.array_equal(added[field],source[field])
                        combined=np.concatenate([original,added]);combined=combined[np.argsort(combined['t'],kind='stable')]
                        assert np.array_equal(modified,combined)
                    elif variant=='outside_context_delete':
                        protected=(abs(original['x'].astype(float)-center[0])<=30)&(abs(original['y'].astype(float)-center[1])<=30)
                        if k>10:
                            previous=old[segment,'F1',k-1];protected|=(abs(original['x'].astype(float)-previous['anchor_x'])<=30)&(abs(original['y'].astype(float)-previous['anchor_y'])<=30)
                        assert np.array_equal(modified,original[protected]) and np.array_equal(removed,original[~protected])
                    else:assert np.array_equal(modified,original)
                    # Independent tuple-set measurement: count occupancy channels, not raw events.
                    def occupied(e):
                        return set((int(x),int(y),bool(p>0)) for x,y,p in zip(e['x'],e['y'],e['p']) if 450<=x<=670 and 120<=y<=280)
                    a,b=occupied(original),occupied(modified)
                    record=next(r for r in interventions if r['segment']==segment and r['variant']==variant and r['frame']==k)
                    assert record['occupancy_changed_polarity_pixels']==len(a.symmetric_difference(b))
                    assert record['occupancy_added_polarity_pixels']==len(b-a)
                    payload_checks.append(dict(segment=segment,variant=variant,frame=k,exact_payload_verified=True,
                                               occupancy_changed_polarity_pixels=len(a.symmetric_difference(b))))
                    if variant=='occlusion' and k in run.ACTIVE:
                        for feature in ['F2','F3']:
                            prior=old[segment,feature,k]
                            if not np.isfinite([prior['anchor_x'],prior['anchor_y']]).all():continue
                            c=[prior['anchor_x'],prior['anchor_y']]
                            n0=int(np.sum((abs(original['x'].astype(float)-c[0])<=13)&(abs(original['y'].astype(float)-c[1])<=13)))
                            n1=int(np.sum((abs(modified['x'].astype(float)-c[0])<=13)&(abs(modified['y'].astype(float)-c[1])<=13)))
                            neighbor.append(dict(segment=segment,frame=k,neighbor_local_feature=feature,original_patch_events=n0,after_target_patch_removal_events=n1,
                                                 retained_fraction=n1/n0 if n0 else np.nan,meaning='patch support only; no physical identity annotation'))
                # Frame13 is the first pressure attempt, frame16 checks removal of the clone itself.
                for k in ([13,16] if variant=='clone_plus8' else [13]):
                    r=next(r for r in rows if r['segment']==segment and r['variant']==variant and r['frame']==k)
                    if 'source_frame' not in r or not np.isfinite(r['source_frame']):continue
                    previous=int(r['source_frame']);source_row=next(x for x in rows if x['segment']==segment and x['variant']==variant and x['frame']==previous)
                    source=data[f'frame_{previous:02d}'];target=data[f'frame_{k:02d}']
                    # Bincount image implementation is independent of B0's scatter accumulation.
                    source_image=run.independent_blurred(source);target_image=run.independent_blurred(target)
                    anchor=np.array([source_row['anchor_x'],source_row['anchor_y']]);template=run.b0.patch(source_image,anchor)
                    scores=np.zeros((25,25))
                    for iy,dy in enumerate(range(-12,13)):
                        for ix,dx in enumerate(range(-12,13)):
                            patch=run.b0.patch(target_image,anchor+[dx,dy]);den=np.linalg.norm(template)*np.linalg.norm(patch)
                            scores[iy,ix]=np.sum(template*patch)/den if den>0 else 0
                    iy,ix=np.unravel_index(np.argmax(scores),scores.shape)
                    assert np.array_equal(anchor+[ix-12,iy-12],[r['proposed_x'],r['proposed_y']])
                    assert abs(scores[iy,ix]-r['score'])<1e-10
                    scores[max(0,iy-3):iy+4,max(0,ix-3):ix+4]=-np.inf
                    assert abs(np.max(scores)-r['second_score'])<1e-10
                    independent_scores.append(dict(segment=segment,variant=variant,frame=k,all_displacements_match=True,second_score_matches=True))
        for variant in run.VARIANTS:
            q=[r for r in rows if r['segment']==segment and r['variant']==variant];lost=False
            for r in q:
                assert r['max_frame_consumed']==r['frame']
                if lost:assert r['status']=='not_attempted_after_loss' and not r['position_valid']
                if str(r['status']).startswith('rejected_'):lost=True
    assert len(contexts)==3*9*4 and all(r['blurred_arrays_bitwise_equal'] for r in contexts)
    return dict(baseline_forward_rows_match_B0_exactly=90,intervention_rows_complete=90,
                previous_B0_files_unchanged=True,run_and_protocol_hashes_match=True,
                intervention_archives_exact_payload_checks=len(payload_checks),occupancy_checks_independent_tuple_sets=True,
                matching_context_access_arrays_bitwise_equal=len(contexts),context_delete_outputs_unchanged=True,
                no_reconnection_all_9_pressure_paths=True,normal_matcher_consumes_current_or_prior_only=True,
                independent_direct_score_maps=len(independent_scores),independent_scores=independent_scores,payload_checks=payload_checks,
                summarize_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),neighbor


def report(baseline,rows,interventions,summaries,validation,neighbors):
    lines=['# B1：布局无模板局部关联压力测试','',
           '结论：B0本来就不需要三维布局模板，但仍依赖局部外观模板。原三段各100ms数据中，本轮只评估参考frame10向未来frame19的50ms。三例目标patch事件移除均在首个干预块拒绝且不复活；+8px同帧复制混合未触发预声明复制锚点吸引，但造成最大1.414px相对基线偏离，其中一例在副本撤去后歧义拒绝。这只是固定输入压力证据，尚不能称可靠拒绝一般干扰、商标识别或跨真实布局泛化。','',
           '没有人工商标/Logo标签，没有真实换布局记录。所有F仍是局部特征名；基线锚点不是物理中心或身份真值。复制相同事件patch也不等同真实新增圆贴点，不能据此否定不规则商标可能提供的额外区分力。','',
           '## 固定协议与评估范围','',
           '预声明见[protocol.md](protocol.md)。保留B0 9种子、5ms、27px双极性binary occupancy外观patch、sigma0.8、±12px整数搜索、12不同像素、score0.35、次峰差0.025、FB2px与不触搜索边的全部门槛。仅前向时间更新外观模板；一次拒绝即停，原信号恢复也不在新位置复用旧F。未增加固定参考模板或亚像素改进。','',
           '9条基线90行的状态、锚点、最佳/次佳分数、FB等匹配字段逐位完全复现B0的5ms前向结果；每F单独重复调用亦相同。它们原本就独立，删除其他seed条目并没有删除邻居事件，也不是换布局实验。前半段仅作原始背景，不拿反向传播补可见性。','',
           '人工参考块为early [1.750,1.755)、middle [6.350,6.355)、late [11.350,11.355)，各向前到1.800/6.400/11.400，共10块含seed、9次自动匹配。正常匹配只用来源和当前事件；FB也只回查来源块，不读取未来块。','',
           '## 输入压力与oracle边界','',
           '固定代表为三个F1，frame13–15即名义中点后15–30ms施加压力。移除/复制位置由B0完整基线轨迹构造，属于非因果oracle输入压力，不用于正常算法初始化/预测/控制。所有原始输入只读，修改后的副本归档为*_events.npz，可逐条核查added/removed/clone_source。','',
           '- `occlusion`文件名仅为简写：实际移除baseline锚点±13px方框内事件，并非遮挡一个已识别物理Marker。方框外原始事件逐条保留，邻居重叠部分可能随方框一起被移除。',
           '- `clone_plus8`复制同帧该方框事件，x+8px，y/t/p及其他字段不改；原事件全保留。两个patch可混合，并非两个独立Marker投影。重复时刻/同像素按B0 binary occupancy归一化，不用raw新增计数冒充有效干扰。',
           '- `outside_context_delete`是有限实现隔离检查：按baseline每帧当前/前一锚点±30px并集保护，删除外部事件。局部外观模板之外，±12px搜索与sigma0.8/truncate4的3px核使源/目标实际上下文至少±28px；30px保守保护并覆盖两相邻步。匹配图像边界x450..670/y120..280不变。此oracle检查不是布局泛化成绩。','',
           '## 干预结果：所有尝试与停止保留','',
           '|片段|输入|有效/10（含seed）|首拒绝帧/原因|接受锚点相对基线最大偏离px|距复制锚点最小值px|停止后帧数|','|---|---|---:|---|---:|---:|---:|']
    for r in summaries:
        rejection='none' if not np.isfinite(r['first_rejected_frame']) else f'f{int(r["first_rejected_frame"])} {r["rejection_reason"]}'
        lines.append(f'|{r["segment"]}|{r["variant"]}|{r["valid_including_manual_seed"]}|{rejection}|{f(r["max_accepted_deviation_from_baseline"])}|{f(r["min_accepted_distance_to_copy_anchor"])}|{r["not_attempted_after_loss"]}|')
    lines+=['','三个patch移除例均在frame13首个受干预块拒绝；frame16原始事件恢复后仍不自动复活。这支持当前三例的拒绝/不接回行为，不能将3/3当部署成功率，亦不能确认单个物理Marker被遮挡。图中邻近结构仍存在，neighbor_patch_counts.csv记录邻居baseline patch计数：并非把整球或整个搜索区抹空。','',
            '复制混合时三个F1在frame13–15仍接受，最大偏离1.414px；早/晚偏离在复制撤去后持续，中段frame16以歧义拒绝。逐帧模板更新可能继承混合造成的偏置，但没有独立身份标签，不能据此判真实identity switch或绝对定位误差。所有接受位置均未满足“距复制锚点≤3px、距基线≥5px且更近至少3px”的预声明copy_anchor_attraction条件；未发生该几何指标不证明一般身份安全。','',
            '以下为首受干预块frame13及副本撤去frame16的实际匹配证据；全部10块、拒绝提议/次峰坐标保存在observations.csv。','',
            '|片段|输入|帧|状态|提议x,y|次佳x,y|score/margin|FB px|','|---|---|---:|---|---|---|---|---:|']
    for r in rows:
        if r['variant']=='outside_context_delete' or r['frame'] not in [13,16] or not np.isfinite(r.get('proposed_x',np.nan)):continue
        lines.append(f'|{r["segment"]}|{r["variant"]}|{int(r["frame"])}|{r["status"]}|{r["proposed_x"]:.0f},{r["proposed_y"]:.0f}|{r["second_x"]:.0f},{r["second_y"]:.0f}|{r["score"]:.3f}/{r["score_margin"]:.4f}|{r["forward_backward_px"]:.3f}|')
    lines+=['','## 原始新增事件与归一化后实际变化','',
            '表中每例为3个活跃干预块的总raw添加/删除事件，以及每块改变的极性×像素occupancy数范围；ON和OFF在同像素分别计通道。这是binary输入实际改变，不是Marker数量。','',
            '|片段|输入|raw新增/删除总计|每块改变occupancy范围|','|---|---|---:|---:|']
    for r in summaries:
        if r['variant']=='outside_context_delete':continue
        lines.append(f'|{r["segment"]}|{r["variant"]}|{r["total_added_events"]:.0f}/{r["total_removed_events"]:.0f}|{r["changed_occupancy_min"]:.0f}–{r["changed_occupancy_max"]:.0f}|')
    lines+=['','`outside_context_delete`三个F1全部输出与基线相同；源/目标正反向所访问的108个blur数组切片逐位相同。其证明这些特定输入路径对实际上下文外变化不敏感，是预期的局部实现隔离属性。不能推导邻近结构（能进入搜索/blur上下文）任意变化也不影响关联，更不能证明真实换布局泛化。','',
            '## 阶段出口与所缺数据','',
            '当前可以保留：无需三维布局模板的人工种子局部关联；三例强patch证据移除能及时拒绝并维持丢失；上下文外更改的有限隔离。不能保留的宣称：已检测商标、已识别物理Marker、真实换布局鲁棒、一般可靠拒绝、三维轴准确。局部外观模板仍存在，逐帧更新产生偏离的机制尚未排除。','',
            '下一步需要人工确认商标/特征身份的连续原始记录、同球或不同球真实改变邻近贴点布局的重复记录，以及短时真实遮挡/出视野/邻居交叉并有人工身份核对的片段。只有这些数据才能区分不规则商标的区分力、真实布局变化与合成叠加。本轮不新增算法，固定参考模板对逐帧更新可作为独立下一轮消融；先保留本轮全部失败与偏离作为对照。','',
            '## 复现与独立检查','',
            f'基线90行、干预90行（3片段×3输入×10帧）完整唯一；90个归档事件payload直接用原SDK绝对时间mask和独立方框定义核对，clone的t/p及其他字段不变，原事件全保留；occupancy差用独立(x,y,polarity)集合核对。{validation["independent_direct_score_maps"]}个关键匹配的25×25位移分数图和次峰以直接内积独立复算，均一致。108个源/目标blur访问切片逐位相同；全部9条压力路径核对拒绝后不重连。旧B0目录文件hash未变，原始输入manifest/part、协议、脚本、副本hash在metadata.json；验证见validation.json。','',
            '```sh','g1_spin/.venv/bin/python g1_spin/experiments/layout_free_local_20261101/run.py','g1_spin/.venv/bin/python g1_spin/experiments/layout_free_local_20261101/summarize.py','```','',
            '图：baseline_review.png原始参考；occlusion_review.png目标patch移除；clone_plus8_review.png合成复制混合。灰+为oracle基线锚点，洋红×为接受锚点，红×为拒绝提议；失去后不补点。表格包含baseline、observations、interventions、summary、context_checks、neighbor_patch_counts、independent_feature_checks。仅写此新目录，无PnP/旋转/训练/联网/硬件，不安装依赖。']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':main()
