# B0：局部特征观测与100ms短时关联

局部观测/短时patch关联已有实证进展：9个人工候选中7条5ms路径覆盖100ms，但2/5ms相对位移差可达9.95px。几何中心、物理身份和共同刚体仍未验证，PnP不启动。这里不是对全部约14个Marker的成功率。固定9个人工patch锚点与门槛的一次运行结果；F1/F2/F3不是MarkerID或几何中心。位置匹配只能表明局部结构在短时内满足本工程gate，不能确认同一物理Marker，更不能证明三维刚体模型。

## 数据与预先固定方法

完整[协议](protocol.md)、[人工种子与依据](seeds.json)、[布局审计](layout_audit.md)。early_turn 1.70–1.80、middle 6.30–6.40、late 11.30–11.40s。2ms各50片、5ms各20片连续不重叠；原始全序列contact sheets在运行前提供父代理审阅，冻结9种子，不替换失败样本。

参考块分别为名义中点之后[+0,+2ms)、[+0,+5ms)，块mid为+1ms/+2.5ms。两时长人工锚点同像素不表示同一瞬时物理点。27×27双极性distinct-pixel occupancy patch，sigma0.8；逐相邻块±12px整数平移搜索，cosine≥0.35、≥12不同像素、最佳与排除邻域3px后的次峰差≥0.025、前后向误差≤2px且不触搜索边。每步更新模板；任一拒绝即该方向停止，不重新接回。

累计模板更新可能漂移或换到相似结构。前后向一致性、次峰分差和非重叠预测支持都只提供自洽性，不是身份真值。原始事件质心不用于驱动匹配；几何模型不用于选点或填轨迹。匹配图像使用x450..670/y120..280的padding，以容纳搜索与patch；原始contact sheets及轨迹审阅主ROI仍为x485..630/y155..250，主ROI外支持可能未显示但在原始匹配计数中保留。

## 每个种子的完整出口

有效数包含1个人工seed；自动匹配数另列CSV。范围为有效原始块的start/end。终止原因不意味着原始结构或物理Marker消失。

|片段|F|T ms|有效/总帧|有效范围s|向前时间终止|向后时间终止|圆/椭圆有效模型数|
|---|---|---:|---:|---|---|---|---:|
|early_turn|F1|2|37/50|1.704–1.778|frame39:rejected_forward_backward|frame1:rejected_ambiguous_peak+forward_backward|10/18|
|early_turn|F1|5|20/20|1.700–1.800|segment_boundary|segment_boundary|3/6|
|early_turn|F2|2|42/50|1.716–1.800|segment_boundary|frame7:rejected_forward_backward|0/2|
|early_turn|F2|5|20/20|1.700–1.800|segment_boundary|segment_boundary|0/1|
|early_turn|F3|2|2/50|1.748–1.752|frame26:rejected_ambiguous_peak|frame23:rejected_ambiguous_peak+forward_backward|0/0|
|early_turn|F3|5|13/20|1.700–1.765|frame13:rejected_ambiguous_peak|segment_boundary|5/8|
|middle|F1|2|29/50|6.300–6.358|frame29:rejected_ambiguous_peak|segment_boundary|27/26|
|middle|F1|5|20/20|6.300–6.400|segment_boundary|segment_boundary|11/16|
|middle|F2|2|34/50|6.332–6.400|segment_boundary|frame15:rejected_ambiguous_peak|12/15|
|middle|F2|5|20/20|6.300–6.400|segment_boundary|segment_boundary|3/3|
|middle|F3|2|7/50|6.344–6.358|frame29:rejected_ambiguous_peak|frame21:rejected_ambiguous_peak|0/0|
|middle|F3|5|19/20|6.300–6.395|frame19:rejected_ambiguous_peak|segment_boundary|1/5|
|late|F1|2|9/50|11.336–11.354|frame27:rejected_forward_backward|frame17:rejected_ambiguous_peak+forward_backward|9/9|
|late|F1|5|20/20|11.300–11.400|segment_boundary|segment_boundary|18/18|
|late|F2|2|16/50|11.322–11.354|frame27:rejected_ambiguous_peak|frame10:rejected_forward_backward|13/14|
|late|F2|5|20/20|11.300–11.400|segment_boundary|segment_boundary|18/19|
|late|F3|2|15/50|11.350–11.380|frame40:rejected_forward_backward|frame24:rejected_forward_backward|2/5|
|late|F3|5|20/20|11.300–11.400|segment_boundary|segment_boundary|16/19|

5ms七条局部结构可覆盖完整100ms：早F1/F2、中F1/F2、晚F1/F2/F3；早F3和中F3有歧义终止。2ms较多提前停止，不能把“2ms匹配未过gate”解释为无Marker信号。晚段5ms仍有可持续局部结构，因此上轮完整球弧失败不等于局部观测完全丢失。

拒绝统计：{"not_attempted_after_loss": 250, "rejected_ambiguous_peak+forward_backward": 3, "accepted": 345, "manual_seed": 18, "rejected_forward_backward": 6, "rejected_ambiguous_peak": 8}。17次拒绝全部涉及歧义分差或前后向门槛，而非本轮score/最小像素/搜索边界门槛。尤其早5ms F3分差0.02469、中5ms F3分差0.02122，均接近预先固定0.025门槛；本轮保留拒绝、不事后放宽。这说明工程gate可能偏严格，不能断言其为物理遮挡。

## 2ms与5ms的位置一致性

只在5ms有效mid取相邻两片均有效的2ms轨迹线性插值，不跨失效段。原始锚点差受seed时差影响；相对位移差只减去最接近名义中点的一个共同有效时刻偏置（逐F固定，不每帧对齐）。以下为中位/P95/最大欧氏差，单位px；共同参考行的零相对差由定义产生，不能算额外验证。

|片段|F|配对数|原始锚点差 中位/P95/最大|固定初始偏置后的位移差 中位/P95/最大|
|---|---|---:|---|---|
|early_turn|F1|14|2.83/6.77/7.49|2.55/6.91/7.67|
|early_turn|F2|17|2.25/8.71/9.55|2.00/8.46/9.30|
|early_turn|F3|0|NA/NA/NA|NA/NA/NA|
|middle|F1|11|1.90/8.13/8.73|1.80/8.42/9.03|
|middle|F2|13|4.03/9.54/9.98|4.10/9.50/9.95|
|middle|F3|2|2.08/2.43/2.47|0.52/0.98/1.03|
|late|F1|4|0.83/1.37/1.46|1.35/1.90/1.95|
|late|F2|6|2.68/3.96/4.03|2.31/4.24/4.56|
|late|F3|6|4.54/6.75/6.79|3.06/5.19/5.25|

## 跨F的重叠诊断

下表只比较同时有效的人工锚点和27×27patch；IoU是正方形patch的几何交并比，不是Marker轮廓重叠率。重叠会共享事件支持，F不能直接当独立Marker约束；即使IoU为零，也不能单凭分散位置证明不同物理身份。

|片段|T|F对|同时有效帧|最小锚点距离px|最大patch IoU|
|---|---:|---|---:|---:|---:|
|early_turn|2|F1/F2|31|23.35|0.08|
|early_turn|2|F1/F3|2|40.05|0.00|
|early_turn|2|F2/F3|2|72.37|0.00|
|early_turn|5|F1/F2|20|23.43|0.08|
|early_turn|5|F1/F3|13|34.06|0.00|
|early_turn|5|F2/F3|13|66.48|0.00|
|middle|2|F1/F2|13|41.00|0.00|
|middle|2|F1/F3|7|41.00|0.00|
|middle|2|F2/F3|7|51.31|0.00|
|middle|5|F1/F2|20|28.44|0.00|
|middle|5|F1/F3|19|37.66|0.00|
|middle|5|F2/F3|19|41.05|0.00|
|late|2|F1/F2|9|45.28|0.00|
|late|2|F1/F3|2|48.70|0.00|
|late|2|F2/F3|2|42.58|0.00|
|late|5|F1/F2|20|43.60|0.00|
|late|5|F1/F3|20|43.68|0.00|
|late|5|F2/F3|20|42.58|0.00|

模型与质心的跨时长比较同存cross_window_positions.csv及diagnostic_distributions.csv；没有可对齐有效中心时不填值。整数位移量化、逐帧模板更新的累积漂移、两种窗口的事件内容差异与局部ROI截断都可能贡献跨窗差；最高约10px差不能当真实Marker中心漂移或中心估计误差。轨迹差不是绝对位置误差，也不能证明共同刚体。2D投影距离可随姿态变化，本轮没有把2D距离恒定作为刚体必要条件。

## 极性、局部中心模型与时间内最大差

质心按不同像素等权计算，ON/OFF/combined分别记录；它们是不同边缘响应的空间位置统计量，不是几何中心或中心估计误差；不能通过平均ON/OFF质心就认定几何中心。圆/椭圆拟合仅对原始不同像素事件位置作候选边缘模型，至少15点、中心bounds锚点±8px、半径/半轴3..18px；至少8/12角bin、最长空角≤90°、median径向残差≤1.5px、RMS≤2.5px，椭圆轴比≤2且不触bounds。无效时正式中心为NaN，原始拟合与拒绝原因仍留在observations.csv。模型角覆盖本身依赖拟合，不能确认物理边界。

下表为有效观测里的ON−OFF质心距离和圆−椭圆中心距离，中位/P95/最大；后一列只比较双方模型gate有效帧。5ms各候选ON/OFF质心距离中位为5.43–16.06px，反映不同极性边缘响应。圆/椭圆同源拟合差小并非独立准确性证据，亦不能排除共同偏差或结构混合。早F2的5ms全部20块可匹配但圆0/20、椭圆1/20有效，尤其说明可匹配patch不等于可用几何中心。

|片段|F|T|ON−OFF n / 中位/P95/最大|圆−椭圆 n / 中位/P95/最大|
|---|---|---:|---|---|
|early_turn|F1|2|37 / 14.41/16.81/17.77|10 / 0.39/0.60/0.61|
|early_turn|F1|5|20 / 13.96/15.79/15.80|3 / 0.10/0.13/0.13|
|early_turn|F2|2|42 / 5.50/7.11/8.04|0 / NA/NA/NA|
|early_turn|F2|5|20 / 5.43/6.58/8.81|0 / NA/NA/NA|
|early_turn|F3|2|2 / 9.87/13.62/14.03|0 / NA/NA/NA|
|early_turn|F3|5|13 / 12.08/16.99/17.21|5 / 0.10/0.18/0.19|
|middle|F1|2|29 / 14.80/16.64/16.73|26 / 0.07/0.18/0.29|
|middle|F1|5|20 / 12.82/16.27/16.31|11 / 0.01/0.09/0.14|
|middle|F2|2|31 / 14.39/18.06/19.54|12 / 0.39/0.88/0.97|
|middle|F2|5|20 / 10.28/13.74/13.97|1 / 0.26/0.26/0.26|
|middle|F3|2|6 / 8.51/14.25/15.65|0 / NA/NA/NA|
|middle|F3|5|19 / 9.42/13.60/13.90|1 / 2.67/2.67/2.67|
|late|F1|2|9 / 14.62/15.31/15.42|9 / 0.19/0.29/0.29|
|late|F1|5|20 / 14.08/16.76/16.96|17 / 0.06/0.20/0.26|
|late|F2|2|16 / 14.27/15.93/16.03|13 / 0.82/1.54/1.55|
|late|F2|5|20 / 14.70/15.57/15.64|18 / 0.52/1.28/1.29|
|late|F3|2|15 / 15.69/17.55/17.63|2 / 0.56/0.93/0.97|
|late|F3|5|20 / 16.06/17.26/17.53|16 / 1.01/1.65/1.73|

diagnostic_distributions.csv还报告各模型/质心相对锚点或彼此的位置差，以及同轨迹这些偏移在时间内的最大两两变化直径offset_change_diameter_px；不是物理中心漂移真值。观测表future_*字段实际指传播方向下一片非重叠块：direction=+1为时间未来，direction=−1为更早；协议中“未来”应按此实现术语澄清。它只在常速预测patch上评分，不搜索、不用于当前接受；加速、形变和遮挡均可降低它，错误持续结构也可升高。

## 阶段出口

- 局部结构可持续追踪：上述七条5ms完整路径仅在人工种子、±12px搜索和当前门槛下成立；两条5ms压力候选及多数2ms路径提前拒绝。不可把缺失补成完整Marker轨迹。
- 几何中心未确认：patch锚点、极性质心、局部圆/椭圆有效模型是不同观测定义；即使模型通过角覆盖gate，也不证明它是球形Marker的投影中心。单独可见弧、重叠事件和未来形变会改变位置。
- 物理身份未验证：F名仅是局部结构标签，累计更新可能跟随边缘而非同一球面点。通过FB与跨窗对齐不证明物理身份。
- 不具备宣称姿态的条件：见[布局审计](layout_audit.md)，已找到10点r20mm模板但未绑定当前约14点球，当前G1相机K亦未找到；本轮不做PnP、三维刚体或角速度结论。

## 复现和独立检查

630观察行唯一且完整；29个代表接受/拒绝匹配用绝对时间boolean选窗、显式像素集合累积及逐位移直接cosine独立核对最佳与次峰，均匹配。36个方向逐条验证丢失后未重连。人工已知(4,3)px位移恢复、相似双峰歧义拒绝、空patch拒绝、重现后不接回均通过。原始210个窗的独立时间检查在preview_metadata.json。完整原始来源/配置/协议/种子/代码hash在metadata.json和validation.json。

```sh
g1_spin/.venv/bin/python g1_spin/experiments/marker_b0_20261101/continuous_preview.py
g1_spin/.venv/bin/python g1_spin/experiments/marker_b0_20261101/run.py
g1_spin/.venv/bin/python g1_spin/experiments/marker_b0_20261101/summarize.py
```

图：continuous_*_{2,5}ms.png为全部连续原始片；trajectories.png为人工锚点路径；actual_frame_review.png为每段5ms第0/5/10/15/19帧的实际事件和当帧有效patch（失效不补点）；cross_window_differences.png为对齐诊断；failure_cases.png按每片段/时长距seed最近拒绝确定性展示，全部17次拒绝和250帧停止后状态在observations.csv。无真实中心/身份标签，无绝对误差、精确率/召回率、三维旋转结论。本轮仅写此新目录，未改旧产物，未安装依赖、训练、联网或操作硬件。
