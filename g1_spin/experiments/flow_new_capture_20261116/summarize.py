"""Compare new-recording initialization routes to existing approximate references."""
import json
from pathlib import Path

import numpy as np

OUT=Path(__file__).resolve().parent
CAPTURE=OUT.parents[2]/'g1_spin/dual_capture/analysis_20261113'


def angle(first, second):
    first=np.asarray(first);second=np.asarray(second)
    dot=np.dot(first,second)/(np.linalg.norm(first)*np.linalg.norm(second))
    return float(np.degrees(np.arccos(np.clip(abs(dot),-1,1))))


def main():
    rows=json.loads((OUT/'results.json').read_text())
    references={r['begin']:r for r in json.loads((CAPTURE/'rgb_period_refined.json').read_text())}
    labels={17:'竖直正转①',19:'竖直正转②',31:'竖直反转',42:'斜轴'}
    details=[]
    for row in rows:
        k=row['time_s']
        ref=references[15 if k<20 else 29 if k==31 else 40]['rpm']
        axis=np.array([0,-1 if k==31 else 1,0] if k<40 else [.65,.76,0],float)
        axis/=np.linalg.norm(axis)
        for name in ['grid','zero','flow','flow_uncorrected']:
            route=row[name]
            if route is None:
                details.append(dict(time_s=k,window_ms=row['window_ms'],method=name,status='initializer_unavailable'))
                continue
            o=route['outcome']
            details.append(dict(time_s=k,window_ms=row['window_ms'],method=name,
                status='accepted' if o['valid'] else o['reason'],rpm=o['rpm'],axis=o['axis'],
                reference_rpm=ref,rpm_difference=o['rpm']-ref,rpm_relative_difference_percent=100*(o['rpm']/ref-1),
                omega_magnitude_rad_s=o['rpm']*2*np.pi/60,
                axis_deviation_deg=angle(o['axis'],axis) if o['axis'] is not None else None,
                score_relative_to_grid=o['contrast']/row['grid']['outcome']['contrast'],
                search_ms=route['search_s']*1000,objective_calls=route['objective_calls'],
                zero_scan_won=route.get('zero_scan_won')))
    lines=['# 新录制光流初始化测试','',
           '2026年11月16日，Tailscale与SSH连接成功后，从Linux取回11月13日RAW的四个既有250ms事件窗口。原始数据不变，文件SHA256与Linux一致。按原球心插值和半径处理，比较严格局部平面、粗空间法向初值、零初值及方向网格。',
           '', '## 原光流前端与粗初值','',
           '| 片段 | 严格平面数 | 严格S可用 | 粗约束数/区域数 | 未补偿粗RPM | 补偿粗RPM | 补偿粗轴 |',
           '|---|---:|---|---|---:|---:|---|']
    for row in [r for r in rows if r['window_ms']==10]:
        c=row['corrected_coarse'];u=row['coarse'];a=c.get('axis')
        axis_text='('+','.join(f'{v:.3f}' for v in a)+')' if a else '无'
        lines.append(f"| {labels[row['time_s']]} | {row['strict']['metadata']['accepted_planes']} | {row['strict']['initializer']['available']} | {c['constraints']}/{c.get('patches',0)} | {u.get('rpm',float('nan')):.1f} | {c.get('rpm',float('nan')):.1f} | {axis_text} |")
    lines+=['','局部邻域使用完整中央10ms的原始整数像素坐标，未先用4000事件下采样。粗约束按拟合时刻的球心与速度构建 r*(s×[nx,ny,0])·omega=vn−n·center_velocity。保留未补偿候选以检查球心修正的影响。所有旋转轴和RGB参考仅在本汇总中用于评价。','',
            '## 各窗口求解结果','',
            '| 片段 | 窗口ms | 方法 | RPM | 对RGB差值 | 对RGB差值% | 对粗参考轴偏差 | 分数/网格 | 状态 |',
            '|---|---:|---|---:|---:|---:|---:|---:|---|']
    for d in details:
        if d['method']=='flow_uncorrected':
            continue
        if 'rpm' not in d:
            lines.append(f"| {labels[d['time_s']]} | {d['window_ms']} | {d['method']} | — | — | — | — | — | {d['status']} |")
            continue
        lines.append(f"| {labels[d['time_s']]} | {d['window_ms']} | {d['method']} | {d['rpm']:.1f} | {d['rpm_difference']:+.1f} | {d['rpm_relative_difference_percent']:+.1f}% | {d['axis_deviation_deg']:.1f}° | {100*d['score_relative_to_grid']:.2f}% | {d['status']} |")
    lines+=['','## 解释范围','',
            '各路线搜索同一窗口、同一组最多4000个事件、同一球面目标与1200RPM上限。flow为沿粗轴扫描幅值后一次自由三维Nelder–Mead，最终轴允许变化。它与此前固定50ms轴估250ms速度的组合不同。strict为原4px/±2ms/0.15std规则，未改门限。',
            '', 'RGB参考来自数秒稳定段周期，轴参考含近似平行像平面的假设，因此表中数值为相对近似参考的偏差。四个既有探索片段未形成新留出验证。多起点G的原生拒绝状态另存，单次路线的共同质量门不能识别所有局部极值。',
            '', '完整轴向量、rad/s、未补偿路线、目标调用数和耗时见summary.json。局部拟合、严格门逐级失败和所有候选见results.json。']
    lines += ['', '## 本次结论', '',
        '原严格前端四窗均为0平面，主要在时间平面残差门退出。粗空间法向方法得到25至35条约束，均覆盖6个区域并达到初始化资格。10/50ms三种搜索接近相同解，短窗速度偏差未被初始化消除。', '',
        '250ms粗光流路线在正转①得到238.27RPM、相对参考轴9.45度，反转得到248.88RPM、7.24度。两者对比分数分别仅为网格高分解的86.03%和93.30%。正转②及斜轴仍收敛到383.35和288.44RPM，与网格相同偏高。初值能够改变局部解，但原评分偏好仍可能指向与参考不符的结果。', '',
        '球心补偿前后初值及最终路线均保存。正转②的未补偿粗轴使幅值扫描退回零向量，最终被拒绝；补偿后进入网格高分解，仍未得到接近RGB的速度。因此球心补偿改善搜索入口不等于改善物理结果。', '',
        '首次试跑的浮点窗口边界导致个别边界事件改变采样，已保留为results_float_boundaries.json。最终采用整数微秒半开区间，12个10/50/250ms网格重跑与历史结果最大相差2.73e-12RPM。独立复算了粗法向全局方程，并核验4个下载文件与Linux SHA256一致。', '',
        '本次只测试四个既有探索窗口。稳定运行、算法接受及较高评分均不等同于轴和速度的精度验证。']
    (OUT/'summary.json').write_text(json.dumps(dict(details=details,
        limitation='Exploratory existing windows, approximate RGB/axis references; geometry and radius not independently calibrated'),indent=2,ensure_ascii=False)+'\n')
    (OUT/'report.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines[:15]))

if __name__=='__main__':
    main()
