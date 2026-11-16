"""Evaluate generated crossings and report known-truth optimizer comparisons."""
from pathlib import Path
import hashlib
import json

import numpy as np
from scipy.spatial.transform import Rotation

OUT = Path(__file__).resolve().parent
rows = json.loads((OUT / 'results.json').read_text())
protocol = json.loads((OUT / 'protocol.json').read_text())
center = np.array(protocol['center'])
radius = protocol['radius']

# This evaluator reconstructs normals from the original stage4 anchor/normal
# definitions without importing the optimizer or current generator.
anchors_xy = np.array([[.6, 0], [0, .6], [-.5, .2], [.2, -.5], [.4, .4], [-.3, -.4]])
normals_xy = np.array([[1., 1.], [1., -.3], [.2, 1.], [.7, -1.], [1., .1], [.4, 1.]])
anchors = np.c_[anchors_xy, -np.sqrt(1 - np.sum(anchors_xy**2, axis=1))]
normals = []
for point, image_normal in zip(anchors, normals_xy):
    jacobian = np.vstack([np.eye(2), -point[:2]/point[2]])
    normal = jacobian @ np.linalg.solve(jacobian.T @ jacobian, image_normal)
    normals.append(normal / np.linalg.norm(normal))
normals = np.array(normals)

checks = []
lines = [
    '# 球面材料边缘合成与优化初值反例',
    '',
    '2026-11-14。复用阶段4的六条有限材料大圆弧，在每个给定转速和窗口长度下重新求整数像素的边缘穿越时刻。采用常轴旋转的解析三角根，并用独立SciPy旋转残差和极性有限差分核验。没有将旧事件时间缩放为新的实录证据。',
    '',
    '半径120px，球心固定，真值轴为(1,2,3)/sqrt(14)。每场景首先以像素规则(x+2y)%5==0留出约20%的事件，再用剩余事件估计。CMAX统一最多5000事件，光流粗候选使用全部训练事件。所有方法输入一致，真值初值单次NM明确作为oracle诊断。',
    '',
    '| 真值RPM | 窗口ms | 转过角度 | 零初值单次NM RPM | 粗光流＋单次NM RPM | 真值初值单次NM RPM | 网格多初值RPM |',
    '|---|---:|---:|---:|---:|---:|---:|',
]
for row in rows:
    candidates = row['candidates']
    speeds = [candidates[name]['metrics']['rpm'] for name in
              ['zero_single_nm', 'coarse_axis_scan_single_nm', 'oracle_true_single_nm', 'grid_multistart']]
    speed_text = [f'{speed:.2f}' for speed in speeds]
    if not candidates['zero_single_nm']['quality']['passed']:
        speed_text[0] += '（拒绝）'
    lines.append(f"| {row['rpm_truth']} | {row['duration_ms']:.0f} | {np.degrees(row['angle_displacement_rad']):.1f}° | " + ' | '.join(speed_text) + ' |')
    events = np.load(OUT / f"{row['scene']}_events.npz")
    truth = np.load(OUT / f"{row['scene']}_truth.npz")
    xy = (events['xy']-center)/radius
    points = np.c_[xy, -np.sqrt(1-np.sum(xy**2, axis=1))]
    material = Rotation.from_rotvec(-(truth['exact_times']-row['duration_ms']/2000)[:,None]*truth['omega']).apply(points)
    scalar = np.sum(material*normals[truth['edge_ids']], axis=1)
    arc_cosine = np.sum(material*anchors[truth['edge_ids']], axis=1)
    rounding = np.abs(events['times_us']*1e-6-truth['exact_times'])
    train, held = events['training_ids'], events['held_ids']
    assert not np.intersect1d(train, held).size
    assert np.isin(events['selected_ids'], train).all()
    assert np.max(np.abs(scalar)) < 1e-12
    assert arc_cosine.min() >= np.cos(.14)-1e-12
    assert rounding.max() <= .500001e-6
    assert np.all(np.diff(events['times_us']) >= 0)
    held_material = Rotation.from_rotvec(-(events['times_us'][held]*1e-6-row['duration_ms']/2000)[:,None]*truth['omega']).apply(points[held])
    true_distance = radius*np.arcsin(np.clip(np.abs(np.sum(held_material*normals[truth['edge_ids'][held]],axis=1)),0,1))
    checks.append(dict(scene=row['scene'],root_residual_max=float(np.abs(scalar).max()),rounding_max_us=float(rounding.max()*1e6),
                       true_held_median_px=float(np.median(true_distance)),true_held_p95_px=float(np.percentile(true_distance,95)),
                       input_sha256=hashlib.sha256((OUT/f"{row['scene']}_events.npz").read_bytes()).hexdigest()))

lines += ['', '| 场景 | 零初值轴误差 | 粗光流初始轴误差 | 粗光流最终轴误差 | 网格最终轴误差 |', '|---|---:|---:|---:|---:|']
for row in rows:
    coarse = np.array(row['coarse']['omega'])
    target = np.array(protocol['axis'])
    initial_angle = np.degrees(np.arccos(np.clip(coarse@target/np.linalg.norm(coarse),-1,1)))
    candidates = row['candidates']
    lines.append(f"| {row['scene']} | {candidates['zero_single_nm']['metrics']['directed_axis_error_deg']:.3f}° | {initial_angle:.3f}° | {candidates['coarse_axis_scan_single_nm']['metrics']['directed_axis_error_deg']:.3f}° | {candidates['grid_multistart']['metrics']['directed_axis_error_deg']:.3f}° |")
lines += [
    '',
    '3000RPM×50ms覆盖2.5圈。零初值NM收敛在6.33RPM，评分相对零运动只提高0.10%，被no_rotation_evidence拒绝。粗光流提供的初轴经幅值扫描和一次NM到3000.05RPM、0.277°，与网格相近。它提供了零初值在较大角位移时会失效、合理方向候选可以帮助跨越局部最优的具体反例。',
    '',
    '1200RPM×50ms覆盖一整圈，四种方法都得到约1202.7RPM和0.36°轴误差。因此一整圈窗口及可见半球变化并非必然导致大偏差。这些因素在真实录像中可能与事件生成、几何误差和纹理结构共同作用，不能仅凭真实250ms失败归因为时间太长。',
    '',
    '六个场景的真值初值NM都会离开精确真值，并获得略高的IWE分数。分数相对真值最多约提高0.194%，最终轴差约0.28–0.78°。整数采样、离散球面网格和Gaussian平滑可使对比度最优点偏离参数真值，本次没有出现真实录像那样数十度的轴偏差。',
    '',
    '| 场景 | 真值在留出事件上的P95边缘平面距离px | 真值初值NM之后的P95距离px | 优化后IWE/真值IWE |',
    '|---|---:|---:|---:|',
]
for row, check in zip(rows, checks):
    oracle = row['candidates']['oracle_true_single_nm']
    lines.append(f"| {row['scene']} | {check['true_held_p95_px']:.5f} | {oracle['metrics']['held_edge_plane_p95_px']:.5f} | {oracle['score_relative_to_truth']:.7f} |")
lines += [
    '',
    '留出距离使用材料edge身份和真值大圆法向，仅用于独立评价，未供给零初值、粗光流或网格估计器。距离只衡量对应edge大圆平面的法向残差，并不包含沿有限弧端点的误差。',
    '',
    '该生成模型是理想几何材料边缘穿越，不模拟DVS对比阈值、多事件尾迹、不应期、噪声、反光、照明或中心运动。六条弧作为独立局部边缘，未为弧端点生成额外边界事件。结果用于优化器和目标函数分析，不代表真实高速球的精度。仅使用一个轴和一套纹理，未测总体失败概率。',
    '',
    '代码及源文件hash见protocol.json，全部候选、拒绝理由、真值指标见results.json，独立事件几何与留出隔离检查见independent_verification.json。',
]
(OUT / 'report.md').write_text('\n'.join(lines)+'\n')
(OUT / 'independent_verification.json').write_text(json.dumps(dict(checks=checks,all_passed=True),indent=2)+'\n')
print('\n'.join(lines[:25]))
