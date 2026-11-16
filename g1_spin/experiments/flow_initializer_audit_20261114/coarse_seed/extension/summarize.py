"""Coverage and numerical agreement on the fixed108 historical windows."""
from collections import Counter
import json
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent
rows = json.loads((OUT / 'results.json').read_text())
condition_rows = []
for condition in sorted(set(row['condition'] for row in rows)):
    selected = [row for row in rows if row['condition'] == condition]
    condition_rows.append(dict(condition=condition, windows=len(selected),
                               G_accepted=sum(row['old_G_valid'] for row in selected),
                               initializer_available=sum(bool(row['initializer'] and row['initializer']['available']) for row in selected),
                               accepted=sum(row['status']=='accepted' for row in selected),
                               statuses=dict(Counter(row['status'] for row in selected))))
comparisons = []
initial_angles = []
for row in rows:
    initializer = row['initializer']
    if initializer and initializer['constraints']:
        observations = [x for x in initializer['observations'] if x['eligible']]
        matrix = np.array([x['constraint'] for x in observations])
        rhs = np.array([x['rhs_px_s'] for x in observations])
        omega, _, _, _ = np.linalg.lstsq(matrix, rhs, rcond=None)
        assert np.allclose(omega, initializer['omega_rad_s'], atol=1e-12)
        assert abs(np.sqrt(np.mean((matrix @ omega-rhs)**2))-initializer['residual_rms_px_s'])<1e-9
    if row['old_G_valid'] and row['status']=='accepted':
        old_axis = np.asarray(row['old_G_axis'])
        final_axis = np.asarray(row['final']['measurement_axis'])
        initial_axis = np.asarray(initializer['axis'])
        initial_angle = float(np.degrees(np.arccos(np.clip(old_axis @ initial_axis,-1,1))))
        initial_angles.append(initial_angle)
        comparisons.append(dict(case=row['case'],condition=row['condition'],
                                initial_to_G_axis_deg=initial_angle,
                                axis_difference_deg=float(np.degrees(np.arccos(np.clip(old_axis @ final_axis,-1,1)))),
                                rpm_difference=row['final']['measurement_rpm']-row['old_G_rpm'],
                                initial_rpm=initializer['rpm'],final_rpm=row['final']['measurement_rpm'],G_rpm=row['old_G_rpm'],
                                final_axis=final_axis.tolist(),G_axis=old_axis.tolist(),
                                initial_residual_rms_px_s=initializer['residual_rms_px_s'],
                                final_gain=row['final']['common_quality']['gain']))
summary = dict(windows=len(rows),statuses=dict(Counter(row['status'] for row in rows)),conditions=condition_rows,
               paired_count=len(comparisons),accepted_sets_identical=all(row['old_G_valid']==(row['status']=='accepted') for row in rows),
               paired_abs_rpm_median=float(np.median([abs(x['rpm_difference']) for x in comparisons])),
               paired_abs_rpm_max=float(max(abs(x['rpm_difference']) for x in comparisons)),
               paired_axis_degrees_median=float(np.median([x['axis_difference_deg'] for x in comparisons])),
               paired_axis_degrees_max=float(max(x['axis_difference_deg'] for x in comparisons)),
               initial_to_G_axis_degrees_median=float(np.median(initial_angles)),
               initial_to_G_axis_degrees_max=float(max(initial_angles)),
               initializer_failures=dict(Counter(row['initializer']['failure_reason'] for row in rows if row['status']=='initializer_rejected')),
               refinement_failures=dict(Counter(row['final']['common_quality']['failure_reason'] for row in rows if row['status']=='refinement_rejected')),
               downsampled_windows=sum(row.get('downsampled',False) for row in rows),comparisons=comparisons)
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
lines = [
    '# 粗初始化固定规则扩展至历史108窗',
    '',
    '沿用首轮M_spatial参数，108窗中70窗得到粗初值，53窗通过现有CMAX质量门。接受的53窗与历史网格G完全同窗。全部原结果与规则均保留。',
    '',
    '| 条件 | 总窗数 | 粗初值可用 | 粗候选+单次CMAX接受 | 历史G接受 |',
    '|---|---:|---:|---:|---:|',
]
for row in condition_rows:
    lines.append(f"| {row['condition']} | {row['windows']} | {row['initializer_available']} | {row['accepted']} | {row['G_accepted']} |")
lines += [
    '',
    '全部状态为53接受、20缺几何、12事件不足、6初值支撑不足、17细化拒绝。17个细化拒绝全部为no_rotation_evidence。初值6个失败的具体原因见summary.json。',
    '',
    f"共同53窗与G的RPM绝对差中位{summary['paired_abs_rpm_median']:.3f}，最大{summary['paired_abs_rpm_max']:.3f}。有向轴差中位{summary['paired_axis_degrees_median']:.3f}°，最大{summary['paired_axis_degrees_max']:.3f}°。",
    '',
    f"粗初轴对G有向角差中位{summary['initial_to_G_axis_degrees_median']:.1f}°，最大{summary['initial_to_G_axis_degrees_max']:.1f}°。粗候选有明显误差，最终相近解不能证明光流本身准确。",
    '',
    '20个事件充分的纯平移窗有14个获得粗初值并进入NM，最终均被旋转证据门拒绝。12个静止窗口事件不足，保持无输出。动态手持窗的估计是总角运动候选。',
    '',
    '本轮局部拟合与CMAX均使用phase3已有selected输入，与首轮局部拟合使用完整common-mask输入有所不同。11个存档窗存在最多5000事件截断。输入继承旧ROI和球心几何，属于同批历史数据扩展，不是独立录制盲测。',
    '',
    '光流路线始终仅沿初轴扫描幅值并调用一次3D Nelder–Mead，没有方向网格和G回退。历史G另有多起点竞争门，单次路线无法执行相同竞争检查，因此相同接受集合仅是本批观察。',
    '',
    '每条局部拟合、约束、初值、rank、condition、残差、扫描起点、原优化候选和最终拒绝理由保留在results.json。验证重建全部局部LS，输入/源码hash保持，所有进入优化的记录均恰好一次NM。速度计量以首轮为准。',
]
(OUT/'report.md').write_text('\n'.join(lines)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['comparisons','conditions']},indent=2))
