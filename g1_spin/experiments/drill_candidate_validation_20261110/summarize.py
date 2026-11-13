"""Descriptive validation with paired axis diagnostics; no model selection loop."""
import json
from pathlib import Path
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
rows=[json.loads(s) for s in (OUT/'results.jsonl').read_text().splitlines()]
assert len(rows)==88 and len({(r['window'],r['method']) for r in rows})==88
reference=json.loads((OUT/'reference.json').read_text());nominal=reference['nominal_same_setting_rpm']

def angle(a,b,unsigned=False):
    a,b=np.asarray(a),np.asarray(b)
    c=float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)))
    return float(np.degrees(np.arccos(np.clip(abs(c) if unsigned else c,-1,1))))

def median(values):return float(np.median(values)) if values else None

def axis_stats(group):
    axes=np.array([r['result']['axis_camera'] for r in group])
    if len(axes)==0:return dict(mean_axis=None,median_deviation_deg=None,max_deviation_deg=None)
    mean=axes.mean(axis=0)
    if np.linalg.norm(mean)<1e-8:return dict(mean_axis=None,median_deviation_deg=None,max_deviation_deg=None)
    mean/=np.linalg.norm(mean);degrees=[angle(a,mean) for a in axes]
    return dict(mean_axis=mean.tolist(),median_deviation_deg=median(degrees),max_deviation_deg=max(degrees),
                abs_z_median=median([abs(a[2]) for a in axes]))

summaries=[];axis_comparison=[]
for run in [1,2,3,4,5]:
    for condition in (['translation','static'] if run==5 else ['rotation']):
        group=[r for r in rows if r['run']==run and r['condition']==condition]
        good_by={m:[r for r in group if r['method']==m and r['result']['valid']] for m in ['baseline','candidate']}
        common=set(r['window'] for r in good_by['baseline'])&set(r['window'] for r in good_by['candidate'])
        for method in ['baseline','candidate']:
            allrows=[r for r in group if r['method']==method];good=good_by[method]
            rpms=[r['result']['spin_rpm'] for r in good]
            devs=[abs(v-nominal) for v in rpms]
            proj=[angle(r['result']['axis_camera'][:2],r['shaft']['direction'],True) for r in good
                  if r['shaft']['valid'] and np.linalg.norm(r['result']['axis_camera'][:2])>.2]
            pair=[r for r in good if r['window'] in common]
            pairproj=[angle(r['result']['axis_camera'][:2],r['shaft']['direction'],True) for r in pair
                      if r['shaft']['valid'] and np.linalg.norm(r['result']['axis_camera'][:2])>.2]
            summary=dict(run=run,condition=condition,method=method,planned=len(allrows),valid=len(good),
                         rpm_median=median(rpms),reference_abs_deviation_mean_rpm=float(np.mean(devs)) if devs and run!=5 else None,
                         reference_scope='same-setting approximation' if run in [1,4] else 'exploratory_unconfirmed_for_this_orientation' if run in [2,3] else 'none',
                         failure_counts=dict(Counter(r['result']['failure_reason'] for r in allrows if not r['result']['valid'])),
                         axis_all=axis_stats(good),axis_common=axis_stats(pair),common_valid=len(pair),
                         shaft_projection_median_deg=median(proj),shaft_projection_count=len(proj),
                         shaft_common_median_deg=median(pairproj),shaft_common_count=len(pairproj),
                         elapsed_median_s=median([r['elapsed_s'] for r in allrows if r['elapsed_s'] is not None]))
            summaries.append(summary)
        for w in sorted(common):
            a=next(r for r in good_by['baseline'] if r['window']==w);b=next(r for r in good_by['candidate'] if r['window']==w)
            axis_comparison.append(dict(window=w,run=run,angle_deg=angle(a['result']['axis_camera'],b['result']['axis_camera'])))
reverse={}
for method in ['baseline','candidate']:
    a=next(r for r in summaries if r['run']==1 and r['method']==method)['axis_all']['mean_axis']
    b=next(r for r in summaries if r['run']==4 and r['method']==method)['axis_all']['mean_axis']
    reverse[method]=angle(a,b) if a is not None and b is not None else None
out=dict(summary=summaries,paired_axes=axis_comparison,run1_run4_directed_mean_axis_angle_deg=reverse)
supplement=json.loads((OUT/'baseline_supplement.json').read_text())
out['standalone_baseline_supplement']=supplement
out['decision']='retain_experimental_candidate_do_not_promote'
(OUT/'summary.json').write_text(json.dumps(out,indent=2,allow_nan=False))

def fmt(v):return '—' if v is None else f'{v:.2f}'
lines=['# 固定候选扩展验证','',
       '44个区间、88条结果；与阶段3分析窗及上一轮窗口对照的事件时间区间不重叠。仍来自相同录制、共享旧轮廓观测，不能称独立采集盲测。', '',
       '**结论：候选保留，不替换生产默认。** 共同有效的run1/4十个窗口中，候选相对同档位参考平均绝对偏差10.47 RPM，基线20.23 RPM；但候选丢失两个原10ms基线可用窗口，且run3杆投影差恶化2.44°，超过预定2°警示线。', '',
       '## 独立基线可用性补充', '',
       '共享20ms几何检查最初使两法共同拒绝10个区间。独立检查发现run1_03565000和run4_00865000的10ms基线本来可用，实际补跑分别输出249.18/237.95RPM；原88条记录不覆盖，补充2条记录单列baseline_supplement.json。',
       '修正后run1基线5/8、候选4/8；run4基线7/8、候选6/8。全部旋转片段两者均22/32，但改善来自run2新增两窗，同时run1/4各丢一窗，不能用总体数量相等掩盖这一点。', '',
       '共同有效轴离散度在四种固定方向均下降；朝相机run2的|z|中位从0.9932到0.9959。run3杆投影差0.23°→2.67°（4个共同有效杆支持窗），run4 1.47°→3.26°（5窗）；前者触发警示，后者也值得保留观察。run3深度倾角没有真值。',
       '纯平移两法均0/8接受非零；静止四窗两法均事件不足，不能称为准确测零。', '',
       '## 共同准备条件下的原始88条对照', '',
       '基线10ms+fixed；候选20ms+median3。参考为另拍同档位246.91RPM，不是逐窗同步真值。run2/3的RPM参考差仅作探索性统计。', '',
       '|run/条件|方法|有效/预定|RPM中位|平均绝对参考偏差RPM|共同有效轴离散中位°|共同杆投影差中位°|',
       '|---|---|---:|---:|---:|---:|---:|']
for r in summaries:
    lines.append(f"|{r['run']}/{r['condition']}|{r['method']}|{r['valid']}/{r['planned']}|{fmt(r['rpm_median'])}|{fmt(r['reference_abs_deviation_mean_rpm'])}|{fmt(r['axis_common']['median_deviation_deg'])}|{fmt(r['shaft_common_median_deg'])}|")
lines += ['',f'run1/run4有向平均轴夹角：基线{fmt(reverse["baseline"])}°，候选{fmt(reverse["candidate"])}°。接近180°支持反转关系，不能替代真实轴标定。', '',
          '轴离散度是相对于各方法自身平均轴的波动，稳定错误轴也能得低值。杆投影为同一事件相机推导的二维诊断，不是三维轴真值；有效杆支持计数见summary.json。',
          '公平性限制：输入为20ms内两种中心共同内区；baseline仅取末10ms，两者分别最多5000事件。几何准备共用20ms覆盖要求，可能拒绝原本10ms可用的基线。补充baseline_eligibility.json专门审计这部分，不把共同准备失败解释为两法本身相同可用率。',
          '离线几何包含未来观测；报告计算耗时不是在线延迟。所有失败保留null，静止事件不足不算测准零。', '',
          '复现：先run.py（拒绝覆盖），再verify.py、summarize.py。输入文件、协议和几何摘要在freeze.json；逐条输入保存在inputs/。']
(OUT/'report.md').write_text('\n'.join(lines)+'\n')
fig,axes=plt.subplots(1,2,figsize=(11,4))
for method,offset,color in [('baseline',-.12,'#2878b5'),('candidate',.12,'#db7b20')]:
    ss=[r for r in summaries if r['run']<5 and r['method']==method]
    axes[0].plot([r['run']+offset for r in ss],[r['rpm_median'] for r in ss],'o-',color=color,label=method)
    axes[1].plot([r['run']+offset for r in ss],[r['axis_common']['median_deviation_deg'] for r in ss],'o-',color=color,label=method)
axes[0].axhline(nominal,color='gray',ls='--');axes[0].set(ylabel='Median accepted RPM',title='Same-setting reference; runs 2/3 exploratory')
axes[1].set(ylabel='Median directed axis deviation (deg)',title='Common valid windows: stability, not accuracy')
for ax in axes:ax.set(xticks=[1,2,3,4],xlabel='Recording');ax.grid(alpha=.2);ax.legend()
fig.tight_layout();fig.savefig(OUT/'comparison.png',dpi=150);plt.close(fig)
print(json.dumps(out,indent=2))
