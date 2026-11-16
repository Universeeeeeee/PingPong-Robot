"""Compare final candidates under the same fixed-geometry event window."""
import json
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent
rows = json.loads((OUT/'results.json').read_text()) + json.loads((OUT/'results_intermediate.json').read_text())
rows.sort(key=lambda row:(row['case'],row['window_ms']))


def angle(a,b):
    a = np.asarray(a)/np.linalg.norm(a)
    b = np.asarray(b)/np.linalg.norm(b)
    return float(np.degrees(np.arccos(np.clip(a@b,-1,1))))


summary = []
for row in rows:
    reference = row['grid']['outcome']
    entry = dict(case=row['case'],window_ms=row['window_ms'],
                 sampled_events=row['used_events'],full_roi_events=row['full_roi_events'],
                 central10ms_events=row['central10ms_events'],
                 initial_rpm=row['initializer'].get('rpm'),initial_axis=row['initializer'].get('axis'),methods={})
    for name in ['grid','zero','flow']:
        result = row[name]
        if result is None:
            entry['methods'][name] = None
            continue
        output = result['outcome']
        item = dict(rpm=output['rpm'],axis=output['axis'],omega_rad_s=output['omega_rad_s'],
                    contrast=output['contrast'],gain=output['gain'],common_valid=output['valid'],
                    native_valid=output.get('native_valid'),native_reason=output.get('native_reason'),
                    contrast_fraction_of_grid=output['contrast']/reference['contrast'],
                    directed_axis_difference_to_grid_deg=angle(output['axis'],reference['axis']),
                    rpm_difference_to_grid=output['rpm']-reference['rpm'],
                    candidate_rotation_degrees=output['rpm']*.006*row['window_ms'],
                    nm_calls=result['nm_calls'],objective_calls=result['objective_calls'],contrast_calls=result['contrast_calls'],
                    zero_scan_won=result.get('zero_scan_won'),selected_seed=result['selected_seed'])
        assert abs(np.linalg.norm(item['axis'])-1)<1e-12
        assert abs(np.linalg.norm(item['omega_rad_s'])*60/(2*np.pi)-item['rpm'])<1e-10
        entry['methods'][name] = item
    summary.append(entry)
(OUT/'summary_all.json').write_text(json.dumps(summary,indent=2)+'\n')
lines = [
    '# 窗口长度、初始化与三维CMAX局部极值',
    '',
    '2026-11-14。四个旧录像代表时刻各比较10/50/100/150/250ms，三种求解均使用完全相同的ROI事件、球面坐标、IWE目标和1200RPM上限。10/50ms三路线近似同解。250ms的零初值单次NM在四窗均停在低分解；粗法向初轴帮助其中两窗达到网格近解，另两窗仍存在明显差异。',
    '',
    '## 输入与假设',
    '',
    '- 中心时刻为run1 1.500s、run2 2.100s、run3 1.500s、run4 2.500s，相对于各RAW第一CD事件。',
    '- 直接读取本地官方解码分片，使用原B1球心和半径，整个250ms保持固定几何，ROI统一0.8r。长窗球心真实变化尚未补偿，本轮用于优化器对照，不能依靠旧转速参考判断物理精度。',
    '- 局部法向初轴仅由中央10ms的完整ROI事件拟合，沿用固定12px/top6patch/五个2ms段/极性分离的M_spatial方案。任何窗口都不输入钻杆或RPM标签。',
    '- 每个10/50/100/150/250ms窗口统一按时间序列最多4000事件。三种求解使用同一组样本。中心时刻、最大转速、NM简单形和终止条件相同。',
    '- G为方向/幅值网格加多起点NM；zero为零向量直接单次NM；flow为粗轴幅值扫描加单次三维NM，最终允许轴变化。它与50ms估轴后固定轴估250ms速度的组合不同。',
    '- 共同质量门保留gain≥0.08、角位移≥0.03rad、收敛、搜索边界与空间退化条件。G原生多起点竞争状态另存。',
    '',
    '## RPM与目标函数',
    '',
    '| 代表窗 | 时长ms | grid RPM | zero RPM | zero对比度/G | flow RPM | flow对比度/G |',
    '|---|---:|---:|---:|---:|---:|---:|',
]
for row in summary:
    g,z,f = (row['methods'][name] for name in ['grid','zero','flow'])
    lines.append(f"| {row['case']} | {row['window_ms']} | {g['rpm']:.2f} | {z['rpm']:.2f} | {z['contrast_fraction_of_grid']*100:.2f}% | {f['rpm']:.2f} | {f['contrast_fraction_of_grid']*100:.2f}% |")
lines += [
    '',
    '## 250ms有向轴与局部极值',
    '',
    '| 代表窗 | grid有向轴 | zero有向轴 | zero对G轴差 | flow有向轴 | flow对G轴差 |',
    '|---|---|---|---:|---|---:|',
]
for row in summary:
    if row['window_ms'] != 250:
        continue
    g,z,f = (row['methods'][name] for name in ['grid','zero','flow'])
    fmt = lambda a:'('+', '.join(f'{v:.3f}' for v in a)+')'
    lines.append(f"| {row['case']} | {fmt(g['axis'])} | {fmt(z['axis'])} | {z['directed_axis_difference_to_grid_deg']:.1f}° | {fmt(f['axis'])} | {f['directed_axis_difference_to_grid_deg']:.1f}° |")
lines += [
    '',
    'run1的flow速度249.27RPM与G250.09RPM接近，但有向轴差27.0°，其对比度为G的97.41%。这给出速度近似、轴仍不同的具体例子。run3粗轴引导到361.11RPM、对G轴差47.46°、对比度为G的78.05%，粗流初值无法保证进入最佳解。',
    '',
    'run2和run4的flow修复零初值局部极值，最终对G轴差均小于0.004°，对比度相对差小于0.0001%。所有250ms flow幅值扫描均选择非零seed，所以上述改进确实使用了粗轴进行搜索。',
    '',
    '60个候选中58个通过当前共同质量门，仅run3的zero100/150ms因no_rotation_evidence被拒绝。20个G原生状态均接受。单次搜索的显著低分解仍可具有正对比度增益且数值收敛，这些门限不能自动识别局部极值。',
    '',
    '## 100/150ms追加对照',
    '',
    '追加完全相同的100/150ms设置。run1/run3/run4的zero最早在100ms观察到明显低分局部解，对应grid候选约164–167°角位移。run2到150ms仍与G同解，在250ms首次出现差异。本轮五种离散窗口不足以确定精确失效阈值。',
    '',
    '100ms和150ms的flow八组全部达到G近解，最大有向轴差0.13°，各对比度/G均接近100%。所有flow幅值扫描都选择非零seed，因此先前短窗试验的零seed绕过现象无法解释这里的改进。250ms的run1/run3仍是明确反例。',
    '',
    '当前共同质量门只判断相对zero增益、角位移、收敛、搜索边界及空间退化。250ms zero四组相对zero图像的增益仍达23.9%、11.3%、24.1%、11.2%，flow run3也有29.8%，均高于8%阈值。它们拥有更高分竞争解这一事实只有多候选比较才能发现。',
    '',
    '## 解释范围',
    '',
    '短窗看到的角位移较小，当前这些样本的局部搜索可从零向量进入网格近解。窗口增大到接近一圈后，三维对比度目标出现对初始化敏感的多个解。粗运动方向可帮助部分样本，但仍需多候选或其他可验证搜索机制。',
    '',
    'G最高IWE也可能有物理模型偏差，因此本轮只证明优化器目标上的局部极值差异。固定中心跨250ms的假设、0.8r ROI、4000采样与旧10ms实验输入不同，不把相对G差值称为真实轴误差。',
    '',
    '完整omega、轴、score、gain、候选角位移、目标调用数、zero_scan_won保存在summary_all.json、results.json和results_intermediate.json。source文件hash和20组样本检查保存在protocol.json、verification.json。',
]
(OUT/'report_all.md').write_text('\n'.join(lines)+'\n')
print('saved20 comparisons; all axis/vector conversions checked')
for duration in [10,50,100,150,250]:
    selected = [r for r in summary if r['window_ms']==duration]
    print(duration,'maxzero_axis',max(r['methods']['zero']['directed_axis_difference_to_grid_deg'] for r in selected),
          'maxflow_axis',max(r['methods']['flow']['directed_axis_difference_to_grid_deg'] for r in selected))
