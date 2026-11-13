# 新录制前的固定评测工具

2026-11-12。只用于离线研究，不连接硬件、不改变g1_spin生产默认。

## 固定方案

- baseline10：10ms窗口中点固定球心。
- smooth20_fallback：20ms球心三点中值平滑/逐事件补偿，几何不覆盖时回退10ms固定球心。
- outer20_fallback：5ms快门、10ms步长提取9个外侧支撑弧几何；可覆盖时20ms平滑补偿，否则回退上一条路线。

CMAX上界3500RPM、内区0.8r、最多5000事件、Gaussian及质量门限与冻结源码一致。回退只在求解前按几何可用性选择，不因输出与参考不同再次求解。三种方案使用各自事件内区，属于整条路线比较；相同路线共享一次计算并显式记录。所有中心仍为离线观测，使用未来支点，不声称在线或因果。

sources.json记录原源码哈希；frozen/保存CMAX、外侧前端、拟合器、回退规则，避免引用历史实验依赖链。freeze.json锁定执行脚本与参数代码；运行前自动校验。未来变更应新建版本，不静默更新此锁。

## 一次评测

从仓库根运行：

```
g1_spin/.venv/bin/python g1_spin/validation/pre_capture_20261112/evaluate.py session/evaluation.json --out session/evaluation_v1
```

输出目录必须不存在。results.json逐窗保存三方案的转速、相机坐标轴、失败原因、回退路线、事件数量及计算时间；summary.json按条件统计覆盖率和失败；provenance.json记录输入哈希、代码锁和环境。保留每条实际执行路线的事件物理ID和中心。耗时区分轮廓前端与CMAX，不含旧几何生成/文件IO，不能当端到端延迟。

输入格式：

- 事件：官方解码manifest.json（parts为结构化NPY字段x/y/t/p，timestamp_unit=us），或小型NPZ字段x/y/t_us/polarity。t/t_us必须是整数微秒，可重复但不可倒序；polarity为0/1。NPZ整份载入，不用于超大文件。manifest路径和geometry_csv相对evaluation.json解析。
- event_origin_us：事件原始时间戳的整数原点；所有窗口和几何时间相对该原点。禁止混用电脑墙钟。
- geometry.csv：t_us,cx_px,cy_px,valid；无效行填valid=0且坐标留空，时间严格递增。不跨30ms几何缺口；缺测不可填零坐标。
- geometry必须已经表达在事件相机像素系与事件时间系。普通相机原图坐标禁止直接填入。
- 每窗明确end_us、radius_px、固定anchor_px与shaft_direction_xy。这一外侧前端使用已知装置粗方向，并不适合无杆方向先验的飞行球直接调用。
- fixed_center_prior仅在独立指定静止中心时使用，不能由算法拒绝自动设置。
- 可选reference_rpm={"kind":"same_setting_approximate","value":246.9113810186355}。kind也可为synchronized_reference/synthetic_truth；输出统一称reference deviation，不认证标注精度。精确轴误差需另有同坐标系轴真值，当前不自动虚构。

参考与上游球心应分别存储：比较三方案时geometry_csv提供冻结的上游几何；外部参考用于检验其误差。如果把外部参考送入算法，应单独命名为reference-assisted实验，不能和事件独立估计混为一谈。

## 同步检查

人工或硬件匹配至少6次独立同步脉冲，填templates/pulses.csv；建议10次以上，覆盖录制开始、中间、结束，用不等间隔序列降低错配。pulse_id必须一一对应，不按算法转速调整时钟。

```
g1_spin/.venv/bin/python g1_spin/validation/pre_capture_20261112/sync.py session/pulses.csv --budget-us 1000 --out session/clock_map.json
```

前2/3脉冲拟合仿射时钟，末1/3至少2个脉冲留出检查，输出漂移ppm和残差。1000us只是明确写出的初始检查预算，不是硬件已达到或算法充分准确的保证。相机若只有60fps帧级标记，时间量化约16.7ms，不能凭插值宣称达到1ms同步。曝光起点/中点、滚动快门和触发传输延迟必须另审。

已经标定投影到事件像素系的参考CSV可转换时间：

```
g1_spin/.venv/bin/python g1_spin/validation/pre_capture_20261112/align_reference.py session/projected_reference.csv --mapping session/clock_map.json --event-origin-us 1000000 --out session/reference_event_time.csv
```

投影参考格式reference_us,cx_event_px,cy_event_px,valid。只允许同步脉冲覆盖区间内映射，禁止外推；未通过同步残差检查拒绝转换。该工具不做跨相机投影，坐标标定必须先独立完成。

## 受控误差实验

```
g1_spin/.venv/bin/python g1_spin/validation/pre_capture_20261112/sensitivity.py --out session/sensitivity_v1
```

固定247RPM、三个轴、两个种子，20ms，14种条件共84次；包括正确几何、0.5/1/2px中心偏移、50/100px/s速度偏差、±1/±5ms参考错位、半径-1/+1/+3px、50%删事件。同一轴/种子除删事件外共享保守内区输入。时间错位作用于具有已知平移/振动的中心轨迹，不能把恒定中心的错位当作同步实验。

这是SciPy独立SO(3)生成的理想表面点轨迹，不是光照/对数阈值事件模拟。结果只判断几何误差敏感性，不等于真实传感器精度或正式验收阈值。丢事件不是参考相机丢帧；参考丢帧由缺口拒绝测试覆盖。

## 验证

```
g1_spin/.venv/bin/python -m pytest -q g1_spin/validation/pre_capture_20261112/test_preflight.py
```

回归用旧数据，不能把回归一致视为新精度证据。实测结果见RESULTS.md；采集与验收见CAPTURE.md。
