"""Summarize exploratory seed tests and verify saved linear systems."""
import json
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent
rows = json.loads((OUT / 'results.json').read_text())


def angle(a, b):
    first = np.asarray(a) / np.linalg.norm(a)
    second = np.asarray(b) / np.linalg.norm(b)
    return float(np.degrees(np.arccos(np.clip(first @ second, -1, 1))))


comparison = []
for row in rows:
    initial = row['initializer']
    final = row['final']
    if initial['constraints']:
        observations = [r for r in initial['observations'] if r['eligible']]
        matrix = np.array([r['constraint'] for r in observations])
        rhs = np.array([r['rhs_px_s'] for r in observations])
        reconstructed, _, rank, singular = np.linalg.lstsq(matrix, rhs, rcond=None)
        assert np.allclose(reconstructed, initial['omega_rad_s'], atol=1e-12)
        residual = np.sqrt(np.mean((matrix @ reconstructed - rhs) ** 2))
        assert abs(residual - initial['residual_rms_px_s']) < 1e-9
    item = dict(case=row['case'], variant=row['variant'],
                initial_rpm=initial.get('rpm'), initial_axis=initial.get('axis'),
                initial_constraint_rms_px_s=initial.get('residual_rms_px_s'),
                final_rpm=final['measurement_rpm'] if final else None,
                final_axis=final['measurement_axis'] if final else None,
                failure_reason=final['common_quality']['failure_reason'] if final else initial['failure_reason'],
                G_rpm=row['old_G_measurement_rpm'], G_axis=row['old_G_axis'],
                pipeline_ms=row['pipeline_s'] * 1000,
                initial_to_G_directed_angle_deg=None, final_to_G_directed_angle_deg=None)
    if row['old_G_axis'] is not None and initial.get('axis') is not None:
        item['initial_to_G_directed_angle_deg'] = angle(initial['axis'], row['old_G_axis'])
    if final and final['measurement_axis'] is not None and row['old_G_axis'] is not None:
        item['final_to_G_directed_angle_deg'] = angle(final['measurement_axis'], row['old_G_axis'])
    if 'G_recomputed' in row and 'search' in row['G_recomputed']:
        g = row['G_recomputed']
        item['G_ms'] = 1000 * g['search']['search_s']
        item['G_recomputed_rpm'] = g['final']['measurement_rpm']
        if row['variant'] == 'clean' and g['final']['measurement_rpm'] is not None:
            assert g['final']['measurement_rpm'] == row['old_G_measurement_rpm']
            assert g['final']['measurement_axis'] == row['old_G_axis']
    comparison.append(item)
(OUT / 'summary.json').write_text(json.dumps(comparison, indent=2) + '\n')
lines = [
    '# 粗法向运动候选能否进入单次CMAX搜索',
    '',
    '2026-11-14。本实验固定一种宽松粗初始化配置，并实际运行光流候选、沿轴幅值扫描和单次三维Nelder–Mead。四个旧旋转窗口均收敛至已有网格法附近的结果，纯平移窗口最终拒绝非零旋转，静止窗口无足够局部支持。',
    '',
    'run1的时间置乱初值仍能被干净CMAX修复至相同解。当前结果支持尝试廉价粗候选搜索入口，无法单独证明粗光流具有可信运动精度或比任意初值更有效。',
    '',
    '## 固定输入与配置',
    '',
    '- 使用phase4存档的完整common-mask事件提取局部运动，保留旧B1球心、半径和区域掩膜。它覆盖旧裁剪区域，不能称为完整球或完整RAW。',
    '- 12px固定网格，按完整10ms事件数量取前6个网格，按2ms分为5段，每个极性单独拟合。局部至少12个独立像素，复用phase5的M_spatial法向速度回归。',
    '- 全局采用r(s×[nx,ny,0])·omega=vn，至少6条观测、3个patch、rank3、condition≤100。球心速度使用B1的0假设。保留全部局部与全局残差，不用此前0.5px最终测量线或残差门限阻止粗候选进入搜索。',
    '- 四个旋转窗及纯平移窗有24–56条观测，rank均3，条件数约2.8–5.0。该条件数只描述已拟局部法向约束，不代表物理可信度。',
    '- 搜索复用phase3 axis_search，3500RPM上限，沿初始有向轴搜索含零的幅值候选，再调用一次3D Nelder–Mead。没有方向网格或G回退。',
    '- CMAX与旧G使用完全相同的phase3 selected事件和球面输入，最多5000事件。run1局部运动输入5460事件，搜索输入5000事件，其余窗口相同。',
    '- 所有轴参考、RPM参考、旧G结果均只用于结果比较，未传入局部拟合和候选搜索。当前六窗口为既有探索样本，没有新增盲测数据。',
    '',
    '## 干净输入的候选与搜索结果',
    '',
    '| 窗口 | 观测数 | 粗候选RPM | 粗候选对G有向轴差 | 全局约束RMS(px/s) | 单次CMAX最终RPM | 网格G RPM | 最终对G轴差 | 候选链/G耗时ms |',
    '|---|---:|---:|---:|---:|---:|---:|---:|---:|',
]
for row, item in zip(rows, comparison):
    if row['variant'] != 'clean':
        continue
    initial = row['initializer']
    if item['initial_rpm'] is None:
        lines.append(f"| {row['case']} | 0 | 无 | — | — | 无 | 无 | — | {item['pipeline_ms']:.1f}/— |")
        continue
    first_angle = f"{item['initial_to_G_directed_angle_deg']:.1f}°" if item['initial_to_G_directed_angle_deg'] is not None else '无G旋转参考'
    final_angle = f"{item['final_to_G_directed_angle_deg']:.3f}°" if item['final_to_G_directed_angle_deg'] is not None else '—'
    final_rpm = f"{item['final_rpm']:.3f}" if item['final_rpm'] is not None else '拒绝'
    g_rpm = f"{item['G_rpm']:.3f}" if item['G_rpm'] is not None else '拒绝'
    lines.append(f"| {row['case']} | {initial['constraints']} | {item['initial_rpm']:.2f} | {first_angle} | {item['initial_constraint_rms_px_s']:.1f} | {final_rpm} | {g_rpm} | {final_angle} | {item['pipeline_ms']:.1f}/{item.get('G_ms',0):.1f} |")
lines += [
    '',
    '以上夹角是对既有G估计的数值差异。G自身不是旋转真值。单次CPU计时包含候选局部拟合及搜索，排除文件读入和原几何提取，属于本次环境测量。四个旋转窗候选链22–36ms，网格G341–506ms，较少的目标函数调用是主要原因。',
    '',
    '纯平移粗候选约63RPM，搜索原候选约58.57RPM，最终因no_rotation_evidence拒绝。不能把粗候选直接当旋转输出。静止仅18事件，全部局部拟合支撑不足。',
    '',
    '## 时间置乱控制',
    '',
    '仅run1和纯平移采用固定seed20261114整体置乱时间，保留xy、极性及时间分布。patch选择仍只由事件计数决定。',
    '',
    '| 窗口与控制 | 粗候选RPM | 最终状态 |',
    '|---|---:|---|',
]
for item in comparison:
    if item['variant'] == 'clean':
        continue
    status = f"{item['final_rpm']:.3f}RPM" if item['final_rpm'] is not None else item['failure_reason']
    lines.append(f"| {item['case']} / {item['variant']} | {item['initial_rpm']:.2f} | {status} |")
lines += [
    '',
    'run1只置乱初始化得到94.47RPM及偏斜轴，干净CMAX仍恢复264.97RPM，与G仅差0.14RPM。完整置乱连同CMAX输入后，单次CMAX和G都拒绝输出。纯平移的两种置乱控制也均拒绝旋转。',
    '',
    '干净四旋转的沿轴扫描均选择350RPM非零起点，因此轴候选确实进入了NM。纯平移和完整时间置乱的扫描选择零起点，此后NM结果不再携带初始轴方向。',
    '',
    '## 结论与证据边界',
    '',
    '旧局部光流严格门限造成没有搜索入口，和粗候选是否足以引导最终CMAX，是两个不同问题。本轮放行粗M_spatial候选后成功执行了真实数据完整链。粗候选本身误差大、全局残差420–1226px/s，最终相似解有明显CMAX纠正贡献。下一步比较应将纯零/随机轴同预算单次搜索纳入，以量化光流是否提供超出廉价任意初始化的收益。',
    '',
    '本次没有独立测定物理轴或精确RPM，没有评估泛化成功率。数学重建检查通过，四个干净G结果与旧存档逐值一致，所有源文件hash前后保持。完整约束、残差、扫描候选、优化器原输出与质量拒绝理由保存在results.json。',
]
(OUT / 'report.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines[25:36]))
