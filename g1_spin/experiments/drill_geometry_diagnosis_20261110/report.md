# 几何缺口与斜轴诊断结果

完成8个既有问题/控制区间 × 2种窗口 × 3种中心，共48条状态；另在前轮44区间上验证只依赖几何的路由。固定范围诊断，不是新独立采集。

## 结论

1. 两个漏窗来自40ms几何观测空洞。20ms不覆盖时回退到有完整几何的10ms固定中心，可以恢复两窗；不填缺口、不放宽30ms限制、不看旋转得分决定回退。
2. 斜轴的补偿几何是优先排查对象。10ms下fixed/linear/median3杆投影差中位约0.24/5.08/2.31度；20ms约1.18/2.27/2.65度。固定时长也能出现差异，不能只归因延长窗口；但无法据此唯一判定真实中心轨迹谁更准。
3. 前后10ms杆投影差约0.25–1.29度；弱参考本身有波动，但这四窗没有显示足以解释所有中心方案差异的统一漂移。杆投影来自同传感器且不是三维轴真值。
4. 不能统一退回固定中心：两个纯平移控制在20ms fixed误报46.85/66.41RPM，而linear/median3未接受非零。需要兼顾平移补偿与轴偏差。
5. 保留实验候选及几何回退，不推广生产默认。后续应优先独立检查中心轨迹/有效观测时刻，不能通过约247RPM择峰或按钻杆方向强行修轴。

## 斜轴六组对照

|窗口ms|中心|有效/预定|RPM中位|杆投影差中位°|
|---:|---|---:|---:|---:|
|10|fixed|4/4|281.84|0.24|
|10|linear|4/4|244.79|5.08|
|10|median3|4/4|264.07|2.31|
|20|fixed|4/4|285.20|1.18|
|20|linear|4/4|270.76|2.27|
|20|median3|4/4|263.85|2.65|

同一时长三种中心共享事件ID；跨时长不同样本且有5000上限，不是纯时长效应。与上一轮差异还含共同掩膜变化，未覆盖旧结果。run3没有独立三维轴和该录制的同步RPM标签。

## 几何路由回放

|run/条件|有效/预定|路由计数|
|---|---:|---|
|1/rotation|5/8|{'geometry_unavailable': 3, '20ms_median3': 4, '10ms_fixed_fallback': 1}|
|2/rotation|8/8|{'20ms_median3': 8}|
|3/rotation|4/8|{'20ms_median3': 4, 'geometry_unavailable': 4}|
|4/rotation|7/8|{'10ms_fixed_fallback': 1, '20ms_median3': 6, 'geometry_unavailable': 1}|
|5/translation|0/8|{'20ms_median3': 8}|
|5/static|0/4|{'20ms_median3': 4}|

旋转24/32可输出（原候选22/32，原基线经补测22/32）；纯平移0/8接受非零；静止0/4因事件不足，不能视作准确测零。所有数字来自记录级复用，不是重新执行完整在线流水线。
回退两个窗口使用前轮独立10ms基线的实际结果249.18/237.95RPM；完整20ms仍使用前轮候选，斜轴问题并未被路由修复。严格检查区间内部缺口，避免只检首尾产生自洽假象。
几何支点含离线未来观测；实时性未验证。路由不会在CMAX不收敛/拒绝之后二次挑结果。

## 核验与复现

5项路由边界测试通过；48状态/14份输入的原始事件ID、公共掩膜、中心插值与平滑独立重建通过；44项路由由逐间隔检查交叉核对；冻结源码及旧实验哈希未变。
`g1_spin/.venv/bin/python -m pytest -q g1_spin/experiments/drill_geometry_diagnosis_20261110/test_route.py`
`g1_spin/.venv/bin/python g1_spin/experiments/drill_geometry_diagnosis_20261110/run.py`（拒绝覆盖）
`g1_spin/.venv/bin/python g1_spin/experiments/drill_geometry_diagnosis_20261110/verify_report.py`
source: protocol.md、freeze.json、geometry_audit.json、routing.json、results.jsonl、verification.json。
