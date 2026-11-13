"""Compare completed fits to approximate references, preserving source hashes."""
import hashlib,json
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parent
files=['spin_10ms.json','spin_short50ms.json','spin_preliminary.json','hybrid_exploratory.json','rgb_period_refined.json']
data={f:json.loads((P/f).read_text()) for f in files}
short={r['time']:r for r in data['spin_short50ms.json'] if 'free' in r}
long={r['time']:r for r in data['spin_preliminary.json'] if 'free' in r}
hybrid={r['time']:r for r in data['hybrid_exploratory.json']}
ten={(r['time'],r['bound_rpm']):r for r in data['spin_10ms.json'] if 'free' in r}
refs={r['begin']:r for r in data['rgb_period_refined.json']}
labels={17:'竖直正转①',19:'竖直正转②',31:'竖直反转',42:'斜轴旋转'}
rows=[];bounds=[]
def unit(a):
 a=np.asarray(a,float);return a/np.linalg.norm(a)
def angle(a,b,unoriented=True):
 dot=float(unit(a)@unit(b));return float(np.degrees(np.arccos(np.clip(abs(dot) if unoriented else dot,-1,1))))
for k in labels:
 ref_rpm=refs[15 if k<20 else 29 if k==31 else 40]['rpm']
 ref_axis=unit([0,-1 if k==31 else 1,0] if k<40 else [.65,.76,0])
 for method,r in [('固定轴250ms',long[k]['one']),('自由三维10ms',ten[k,1200]['free']),('自由三维50ms',short[k]['free']),('自由三维250ms',long[k]['free']),('组合50+250ms',hybrid[k]['result'])]:
  if 'signed_rpm' in r:
   rpm=abs(r['signed_rpm']);axis=unit(r['axis'])*(1 if r['signed_rpm']>=0 else -1);valid=True
  else:
   valid=r['valid'];rpm=r.get('spin_rpm');axis=r.get('axis_camera')
  row=dict(time_s=k,phase=labels[k],method=method,valid=valid,reference_rpm=ref_rpm,reference_axis=ref_axis.tolist(),reference_axis_kind='manual coarse image-plane assumption')
  if rpm is not None and axis is not None:
   axis=unit(axis);speed=rpm*2*np.pi/60;ref_speed=ref_rpm*2*np.pi/60
   row.update(axis=axis.tolist(),axis_deviation_deg=angle(axis,ref_axis),directed_axis_deviation_deg=angle(axis,ref_axis,False),rpm=rpm,omega_magnitude_rad_s=speed,omega_vector_rad_s=(axis*speed).tolist(),rpm_difference=rpm-ref_rpm,rpm_relative_difference_percent=100*(rpm/ref_rpm-1),omega_magnitude_difference_rad_s=speed-ref_speed)
  else:row['failure_reason']=r.get('failure_reason')
  rows.append(row)
 a=ten[k,1200]['free'];b=ten[k,3500]['free']
 bounds.append(dict(time_s=k,rpm_1200=a.get('spin_rpm'),rpm_3500=b.get('spin_rpm'),axis_difference_deg=angle(a['axis_camera'],b['axis_camera']) if a['valid'] and b['valid'] else None))
summary=dict(rows=rows,bound_check=bounds,excluded=[dict(time_s=45,reason='insufficient_geometry; only 2 unambiguous circle detections')],sources={f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in files},limitations=['Approximate segment RGB speed reference, not instantaneous encoder ground truth','Manual coarse shaft axes with image-plane assumption, not calibrated 3D truth','Four exploratory windows; no held-out validation','10/50/250ms centered at k+0.125s; common geometry and 0.8r mask, deterministic up to 4000 events; changing span changes actual event sample','10ms 3500RPM bound matches earlier bound but does not match old recording/geometry/event preprocessing'])
(P/'methods_comparison.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
lines=['# 本次录制：10ms补跑与方法汇总','', '10ms自由三维补跑完成。主表统一1200RPM上限；额外3500RPM上限检查用于对照上周搜索范围。轴为有向单位向量，x右、y下、z远离相机。角速度大小单位rad/s。RPM差值为估计减RGB段参考；角度偏差按无向轴线夹角计算。固定轴方向由人工提供，其0°为定义结果。','', '| 片段 | 方法 | 轴(x,y,z) | 相对参考轴偏差 | RPM | 角速度大小rad/s | RPM差值 | RPM相对差值 |','|---|---|---|---:|---:|---:|---:|---:|']
for r in rows:
 a='('+', '.join(f'{x:.3f}' for x in r['axis'])+')' if 'axis' in r else '拒绝输出'
 lines.append(f"| {r['phase']} | {r['method']} | {a} | {r.get('axis_deviation_deg',float('nan')):.1f}° | {r.get('rpm',float('nan')):.1f} | {r.get('omega_magnitude_rad_s',float('nan')):.2f} | {r.get('rpm_difference',float('nan')):+.1f} | {r.get('rpm_relative_difference_percent',float('nan')):+.1f}% |")
lines+=['','## RGB近似参考','', '| 片段 | RGB参考RPM | 对应角速度rad/s |','|---|---:|---:|']
for k in labels:
 r=next(r for r in rows if r['time_s']==k);lines.append(f"| {labels[k]} | {r['reference_rpm']:.2f} | {r['reference_rpm']*2*np.pi/60:.2f} |")
lines+=['','## 搜索上限检查','', '| 片段 | 1200上限RPM | 3500上限RPM | 两次估轴夹角 |','|---|---:|---:|---:|']
for r in bounds:lines.append(f"| {labels[r['time_s']]} | {r['rpm_1200']:.3f} | {r['rpm_3500']:.3f} | {r['axis_difference_deg']:.3f}° |")
lines+=['','## 解释与范围','', '10ms轴方向总体接近粗参考，但RPM相对差值约−21%至+16%；50ms并非每个片段都更好。250ms自由三维转速与轴均出现较大偏离。组合在这四个窗口的速度参考差值最小，但它复用了50ms轴，轴精度并没有再次改善。','', '1200与3500上限的10ms结果最大差值小于0.1RPM，未解释本批速度偏差。这只匹配上周窗口长度和搜索上限，不能视作完全匹配上周全部处理设置。','', '45s预选窗口因几何检测不足而排除。所有误差均是相对近似参考的偏差：RGB是数秒稳定段内周期估计，轴参考包含平行像平面假设。四窗口为探索结果，不能外推为真实击打精度或整体成功率。','', '本次10ms补跑使用的时间区间为k+0.120至k+0.130秒，与已有50ms、250ms共用中心时刻k+0.125秒。所有方法保持已有球心插值、半径、0.8r球面区域及最多4000事件采样规则。','', '完整数值（含角速度向量、rad/s差值和源文件SHA256）见methods_comparison.json。']
(P/'methods_comparison.md').write_text('\n'.join(lines)+'\n')
assert len(rows)==20 and all(r['valid'] for r in rows)
assert all(abs(np.linalg.norm(r['axis'])-1)<1e-10 for r in rows)
assert all(abs(np.linalg.norm(r['omega_vector_rad_s'])-r['omega_magnitude_rad_s'])<1e-10 for r in rows)
assert all(hashlib.sha256((P/f).read_bytes()).hexdigest()==h for f,h in summary['sources'].items())
print('\n'.join(lines[:25]));print('VERIFIED: 20 rows; unit axes; angular-speed vectors consistent; original sources unchanged')
