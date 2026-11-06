"""Summarize fixed protocol results and independently check crop/time bookkeeping."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates
import run

OUT=Path(__file__).resolve().parent


def read(name):
    return list(csv.DictReader((OUT/name).open()))


def values(rows,key):
    return np.array([float(r[key]) for r in rows])


def main():
    rows=read('observations.csv')
    paired=read('paired_windows.csv')
    pert=read('block_perturbations.csv')
    radii=read('radius_sensitivity.csv')
    summary=[]
    for segment in run.SEGMENTS:
        for duration in [20,50,100]:
            for mode in ['free','fixed48']:
                rr=[r for r in rows if r['segment']==segment and r['duration_ms']==str(duration) and
                    r['mode']==mode and r['threshold_policy']=='scaled_0.6']
                ss=dict(segment=segment,duration_ms=duration,mode=mode,n=len(rr))
                for key in ['cx','cy','radius','rms','median_abs_residual','candidate_count',
                            'candidate_raw_2px_isolated_fraction','supported_2px_longest_original_rayrun_deg',
                            'occupied_15deg_bins','supported_2px_15deg_bins','raw_circle_2px_15deg_bins',
                            'raw_lower_semicircle_2px_bins','jac_condition','jac_cy_r_cos']:
                    ss['median_'+key]=np.median(values(rr,key))
                for key in ['cx','cy','radius']:
                    ss['range_'+key]=np.ptp(values(rr,key))
                    ss['adjacent_difference_rms_'+key]=np.sqrt(np.mean(np.diff(values(rr,key))**2))
                summary.append(ss)
    run.write_csv('summary.csv',summary)
    # Independent full-sensor accumulation cross-checks the cropped implementation.
    manifest=json.loads((run.DATA/'manifest.json').read_text())
    events=np.concatenate([np.load(run.DATA/f) for f in manifest['parts']])
    t0=int(events['t'][0])
    old=np.genfromtxt(run.ROOT/'g1_spin/outputs/real_222403/qc/circle_fits.csv',delimiter=',',names=True)
    checks=[]
    for segment,span in run.SEGMENTS.items():
        mid=sum(span)/2
        for duration in [20,50,100]:
            lo=t0+round((mid-duration/2000)*1e6)
            hi=t0+round((mid+duration/2000)*1e6)
            selected=events[(events['t']>=lo)&(events['t']<hi)]
            c=np.bincount(selected['y'].astype(int)*manifest['width']+selected['x'],
                          minlength=manifest['height']*manifest['width']).reshape(manifest['height'],manifest['width'])
            sm=gaussian_filter(c.astype(float),1)
            profiles=map_coordinates(sm,[run.YY,run.XX],order=1)
            pts=[]
            for i,profile in enumerate(profiles):
                ix=np.flatnonzero(profile>=.6*duration/100)
                if len(ix):pts.append([run.XX[i,ix[-1]],run.YY[i,ix[-1]]])
            points=np.array(pts)
            cropped,_,_=run.candidates(selected,duration,'scaled_0.6')
            row=next(r for r in rows if r['segment']==segment and r['duration_ms']==str(duration) and
                     r['threshold_policy']=='scaled_0.6' and r['mode']=='free' and abs(float(r['midpoint_s'])-mid)<1e-7)
            p,_=run.fit_circle(points)
            check=dict(segment=segment,duration_ms=duration,events_direct_mask=len(selected),
                       event_count_matches=len(selected)==int(row['events_total']),
                       max_fullframe_crop_point_difference=float(np.max(abs(points-cropped))),
                       max_recomputed_fit_csv_difference=float(np.max(abs(p-np.array([float(row[k]) for k in ['cx','cy','radius']])))))
            if duration==100:
                previous=old[np.argmin(abs(old['relative_s']-(mid-.05)))]
                check['old_start_s']=float(previous['relative_s'])
                check['old_pipeline_max_cxcyR_difference']=float(np.max(abs(p-np.array([previous[k] for k in ['cx','cy','radius']]))))
            checks.append(check)
    validation=dict(expected_grid_rows=4*31*3*4*2,actual_rows=len(rows),
                    unique_grid_rows=len({tuple(r[k] for k in ['segment','midpoint_s','duration_ms','threshold_policy','mode']) for r in rows}),
                    optimizer_valid_rows=sum(r['optimizer_valid']=='True' for r in rows),
                    radial_boundary_hit_rows=sum(int(r['radial_boundary_hits'])>0 for r in rows),
                    checks=checks,summarize_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    run_hash_matches_metadata=hashlib.sha256((OUT/'run.py').read_bytes()).hexdigest()==json.loads((OUT/'metadata.json').read_text())['script_sha256'])
    assert validation['actual_rows']==validation['expected_grid_rows']==validation['unique_grid_rows']
    assert validation['run_hash_matches_metadata']
    assert all(c['event_count_matches'] and c['max_fullframe_crop_point_difference']==0 and
               c['max_recomputed_fit_csv_difference']<1e-5 for c in checks)
    assert all(c.get('old_pipeline_max_cxcyR_difference',0)<.001 for c in checks)
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    lines=['# DVSLume-G1 几何观测离线实验（阶段 A）','',
           '结论：现有上弧外侧采点只能产生待筛选的几何候选，阶段 A 尚不足以作为经过验证的瞬时球心观测出口。早/中段50–100ms候选较连续，但累积包络与瞬时轮廓之间仍无独立真值；晚段应退出当前无条件几何观测路径。最强可验证发现是20ms主阈值容许单事件噪声通过，问题属于本轮“旧外侧取点+该阈值”流程，不能由此断言G1本身不能20ms定位。50ms/fixed48也只是候选，未宣布通过。此轮检验旧外侧提取器在各窗下的观测支持，未给出G1硬件的时间分辨率或精度极限，未对轮廓来源作确定分类。没有绝对准确率、绝对像素误差或旋转结论。','',
           '## 固定协议与假设','',
           '- 数据为官方SDK解码的6,020,744条CD事件，17.579495s。相对第一事件197754251us，保持原时间与顺序，未修改输入；manifest及25个part的SHA256见metadata.json。',
           '- 片段：平移对照0.30–0.60s、早段转折1.60–1.90s、中段6.20–6.50s、晚段11.20–11.50s。每段31个统一名义中点，步长10ms；20/50/100ms窗均为[mid−T/2,mid+T/2)，从片段边界外取足数据。start/end与mid均保留。',
           '- 原点(556,224)、上半圆121射线、r=30..75.5步长0.5、sigma=1 blur、取每射线最后一个过阈值点均保留旧流程；原点只负责采样/初始化，不是真值。soft-L1径向拟合，原有中心/半径bounds不变。主表free与固定r48；代表中点另试r46/50。',
           '- 主阈值为0.6×T/100ms，保持平均事件密度阈值；对照固定0.6 count，以及0.3/1.0×T/100ms。阈值协议在查看结果前固定。20ms主阈值0.12低于单个事件的高斯峰约0.159，故单个外噪事件可触发；这不是等物理信息量比较。',
           '- 原始支持ROI包含下半圆，采样ROI仅包含上弧；未用下方截断边拟合。lower_occlusion_risk是全片既有遮挡先验标记，不是本脚本自动识别的结论。',
           '- “有效角度”同时记录候选15°bin占用、拟合残差≤2px的bin占用/最长相邻bin run，及更严格的连续原射线run。后者单位为原点角度，不能与拟合球心角度直接互换。bin内一个点不能证明整段连续；采到121射线也不证明真实轮廓。','',
           '## 主阈值的分段结果','',
           '下表均为free fit，31窗中位数；“孤点”为候选点2px邻域内原始事件≤1的比例。RMS包含离群候选，soft-L1目标并非最小化RMS。','',
           '|片段|窗ms|cx范围|cy范围|r范围|RMS中位px|孤点比例|连续射线run中位°|','|---|---:|---:|---:|---:|---:|---:|---:|']
    for s in summary:
        if s['mode']=='free':
            lines.append(f"|{s['segment']}|{s['duration_ms']}|{s['range_cx']:.2f}|{s['range_cy']:.2f}|{s['range_radius']:.2f}|{s['median_rms']:.2f}|{s['median_candidate_raw_2px_isolated_fraction']:.1%}|{s['median_supported_2px_longest_original_rayrun_deg']:.1f}|")
    lines+=['','cx/cy/r范围含真实平移、采样变化和污染，不是观测误差或噪声标准差。summary.csv另存相邻10ms差分RMS，只是混合变化量：相邻窗重叠50%/80%/90%，不能当独立重复或据此证明长窗更准确。','',
            '## 统一中点的窗口配对与阈值混淆','',
            '|片段|模式|20−100ms绝对dx中位/P95|绝对dy中位/P95|','|---|---|---:|---:|']
    for segment in run.SEGMENTS:
        for mode in ['free','fixed48']:
            q=[r for r in paired if r['segment']==segment and r['mode']==mode and r['short_ms']=='20' and r['long_ms']=='100']
            dx=abs(values(q,'dx_short_minus_long'));dy=abs(values(q,'dy_short_minus_long'))
            lines.append(f'|{segment}|{mode}|{np.median(dx):.2f}/{np.quantile(dx,.95):.2f}|{np.median(dy):.2f}/{np.quantile(dy,.95):.2f}|')
    lines+=['','20ms free fit从主阈值0.12改固定0.6时，四段cy配对绝对差中位依次为'+
            ' / '.join(f"{np.median(abs(values([r for r in rows if r['segment']==segment and r['duration_ms']=='20' and r['mode']=='free' and r['threshold_policy']=='scaled_0.6'],'cy')-values([r for r in rows if r['segment']==segment and r['duration_ms']=='20' and r['mode']=='free' and r['threshold_policy']=='fixed_0.6'],'cy'))):.2f}" for segment in run.SEGMENTS)+
            'px。时长与阈值/稀疏支持相互影响，不能把窗口差全部归因为运动模糊。100ms主/固定阈值本来相同。','',
            '## cy/r耦合、区域扰动与数值收敛','',
            '- 自由模型Jacobian按[-ux,-uy,-1]构造，jac_condition使用实际参数列；fixed48只用前两列。released_radius_jac_condition与jac_cy_r_cos在固定模型中表示“若释放半径”的潜在几何耦合，不是固定模型条件数或统计协方差。',
            '- 自由模型各段Jacobian列cy/r余弦中位约−0.84至−0.89，表明上弧的纵向移动与半径变化部分互相补偿；三初值收敛近一致不能消除几何耦合，也不能证实测量准确。',
            '- 代表中点r46→50时，cy跨度：'+
            '; '.join(segment+' '+','.join(f"{np.ptp(values([r for r in radii if r['segment']==segment and r['duration_ms']==str(d)],'cy')):.2f}" for d in [20,50,100])+'px' for segment in run.SEGMENTS)+
            '（依次20/50/100ms）。固定r48的更稳cy部分来自先验，不能作为精度提升证据。',
            '- 代表中点做四种删除一个连续时间四分块的重拟合（阈值按剩余75%曝光缩放），以及删除一个连续45°原射线区域的重拟合。这些是依赖结构下的敏感性诊断，不是独立bootstrap或置信区间。','',
            '|片段|free时间块最大中心偏移|free空间区最大中心偏移|fixed48空间区最大中心偏移|','|---|---:|---:|---:|']
    for segment in run.SEGMENTS:
        def maxshift(mode,kind):
            q=[r for r in pert if r['segment']==segment and r['mode']==mode and r['perturbation']==kind]
            return max(np.hypot(values(q,'dx'),values(q,'dy')))
        lines.append(f"|{segment}|{maxshift('free','delete_time_quarter'):.2f}|{maxshift('free','delete_45deg_sector'):.2f}|{maxshift('fixed48','delete_45deg_sector'):.2f}|")
    lines+=['',f"主格点完整保留{len(rows)}行（4×31×3×4×2），optimizer_valid={validation['optimizer_valid_rows']}仅指收敛且未触参数bounds，不能称轮廓可靠率。4行无效均为20ms固定0.6阈值的free模型触界（中段6.23s，晚段11.34/11.35/11.36s）。径向采样触顶/触底另有{validation['radial_boundary_hit_rows']}行，详见radial_boundary_hits。",'',
            '## 原始证据与阶段A分段出口','',
            '请先看raw_events_only.png，再看raw_candidate_overlays.png。图中每条原始事件均绘制，未随机抽样；重复像素重叠。ON橙/OFF蓝，候选红圈，自由圆青/固定48绿，黑叉为采样原点；灰线y=240只是既有遮挡关注位置。','',
            '|片段|视觉证据等级与出口判断|','|---|---|',
            '|平移对照|100ms形成更连续的外包络，但平移和纹理扫掠同时累积；可作运动敏感性对照，不能出口为瞬时球心。20ms有外层孤点。|',
            '|早段转折|50–100ms上弧候选较连续，适合作人工几何ROI候选；20ms孤点与固定阈值敏感性明显。尚不能判定silhouette身份。|',
            '|中段|50–100ms有可用上弧候选；20ms free cy/r跨度明显放大，不能无筛选使用。长窗平滑不能单独构成通过证据。|',
            '|晚段|候选明显混入内部纹理与外层孤点，部分弧缺失；固定半径也未消除阈值/区域敏感性。当前阶段A不通过，应退出无条件几何观测路径。|','',
            '下方可见事件在约y=240附近截断且拟合圆底部多数无支撑，与遮挡先验相符；纹理事件也可在圆周附近形成高计数。原始±2px圆周角bin有事件并不能确认轮廓；孤点比例低也可能来自重复噪声/纹理。此处视觉等级是对有限代表图的人工审查，未自动分类全部窗。','',
            '未解决：如何区分真实silhouette、内部纹理扫掠和孤点；如何校准半径先验及cy/r耦合；如何量化时间窗内运动偏差；如何在无独立真值条件下定义可接受的几何观测支持。建议先在早/中50–100ms候选上人工筛选和设计独立几何验证，暂不推进旋转估计。','',
            '## 复现与验证','',
            '在项目根目录运行：','',
            '```sh','g1_spin/.venv/bin/python g1_spin/experiments/geometry_observation_20261101/run.py',
            'g1_spin/.venv/bin/python g1_spin/experiments/geometry_observation_20261101/summarize.py','```','',
            '使用已有numpy/scipy/matplotlib，不安装依赖。仅离线读取解码数据并写本目录；未改s-CMax、optical_flow或旋转估计，无远程连接、训练或硬件控制。','',
            '验证：12个代表窗用全1280×720帧重新累积，独立绝对时间boolean mask核对取窗；与裁剪采样所得候选点逐点完全一致，重拟合与CSV差<1e−5px。4个100ms代表窗与旧流程对应start=mid−50ms的cx/cy/r差<0.001px（见validation.json），支持旧外侧基准被保留。上述仅验证实现一致性，不是几何真值验证。','',
            '文件：observations.csv完整主表；summary.csv分组中位/范围与相邻差分；paired_windows.csv配对时长差；radius_sensitivity.csv参考半径对照；block_perturbations.csv时间/空间删除；multistart.csv三初值；metadata.json输入与协议；validation.json实现检查。图像另有center_radius_traces.png与threshold_cy_sensitivity.png。']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:v for k,v in validation.items() if k!='checks'},indent=2))


if __name__=='__main__':
    main()
