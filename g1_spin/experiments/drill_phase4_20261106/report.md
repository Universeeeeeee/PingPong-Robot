# 阶段4最终报告：合成接口贯通，语义边界明确，实录前端仍未获得可验角速度

**三个任务已按固定范围完成并停止。** task1补齐了同一整数事件贯通local motion、满秩初始化和CMax搜索的正例；task3澄清A_literal解释、平移约束及加权LS边界；task2完成三种前端的严格T/H/V开发、确认、27组合split诊断及108窗回归，仍没有得到可用实录A/S初始化或可验角速度。没有因结果不好更改门限、扩窗或增加第4前端。

本阶段由GPT-6 Astra xhigh审阅、GPT-6.1 Sol high执行，最终交Astra终审。生产源码、阶段1/2/3结果、RAW和task1+3封存文件均保持hash一致。全部实录几何和旋转仍属hypothesis_only、physics_validated=false；同一五段历史录制的回归不能当独立测试。

## task1：同一events全链正例已经通过

独立SciPy Rotation和整数像素过局部材料大圆弧的根求解生成6组固定控制，标签和omega只交evaluator，runner仅接events+geometry。三旋转正例的S完整链均满足预定合成门槛，反向符号正确：

|正例|S初轴差 / speed差|S最终轴差 / speed差|A_literal初轴差→最终轴差|
|---|---:|---:|---:|
|base (10,20,30)|0.067° / +0.110%|0.682° / +0.668%|3.817°→0.710°|
|reverse|0.123° / −0.130%|0.719° / +0.697%|3.336°→0.728°|
|second (−20,30,10)|0.024° / −0.052%|0.120° / −0.028%|6.964°→0.104°|

每个A/S真实调用一次NM，没有trueinit或G回退。S三个初值均比细化后更接近合成参数，尽管IWE提高；不能将目标分数改善泛称精度提高。独立root/极性/球面/LS/IWE核验、未拟合事件预测和跨物理edge留一评价通过，详见[task1+3检查点报告](report_task13.md)、[验证](verification_task13.json)及[合成评价](synthetic_evaluation.json)。它是干净几何edge-crossing实现测试，不含完整DVS噪声/对比阈值/不应期模型，不等于实录成功或论文精度复现。

## task3：A_literal保留，不能臆测官方未公开语义

[semantics.md](semantics.md)记录一次有限来源核查。主文公开plane-fit二维flow→半径归一化→切向lift→P→LS流程，却未明确g到三维normal velocity的细式；官方公开approximate pseudo code有时间表面/CNN/CMax，没有offline OF/P推导；可得补充不足以消歧。不能把官方在线时间表面当offline规定，更不能宣称Sony未公开实现错误。

阶段3原解释保留名A_literal及非等价反例。raw法向约束为`r(s×[gx,gy,0])·omega=1−gᵀvc`；整个events坐标先减中心轨迹后重新提取，是另一明确路径，禁止用`f_normal−vc`冒充完整flow。切向协向量h及其最小范数速度的逐条约束对应S，但stack投影LS隐含1/||h||权重，通常不等于未加权S；b_rhs=0时仍保留h/P，不除零。解析语义测试通过，没有静默替换A_literal或增加第4实录算法。

## task2：前端、预测和选择均在实录前冻结

三前端仅F0全事件4px/±2ms、F1同极性因果latest surface、F2在F1上最多64个seed锚定triplet的一次稳健簇/一次LS重拟合；最终沿用原邻居、空间rank、plane RMS和初始化门限。使用完整stage3 common_ids，先按前6ms及固定2px棋盘分T/H，后4ms为V；所有seed、邻域、簇、LS、CMax只读T。

实际规则以[协议顶部审阅修订](real_frontend/protocol.md)及[parameter_freeze](real_frontend/parameter_freeze.json)为准：先在全部T因果history取每sensor pixel latest，再查lane半径；球面非法source明确失败；背面预测保留3px惩罚；两个方向分别判距≤1.5px/改善≥20%，等权报告双向mean Chamfer；非零旋转必须同时优于zero和translation-only；H还要求≥20独立事件、总H≥20%和三区各≥5事件；被拒绝LOO候选不计通过。六个实际路径控制均通过，初版仅合成输出/代码快照保留，没有实录后调参。

同样六组合成scene在新T/H/V下，3旋转×3F的9个S正例均通过全链与双方向V线；这说明F1/F2实现可处理干净模型。慢纯平移的translation-only却未通过V工程线（两向中位1.793/1.587px，forward改善19.81%），揭示固定1ms模板/离散稀疏采样的指标局限，不能据此证伪已知运动，也没有为该控制放宽门槛。

### 开发、确认与split测量诊断

4开发窗（run1 1.495/run2 2.095/run5 4.095/静止1.525）和4确认窗（run1 1.425/run3 1.495/run4 2.495/run5 5.095）均保留3F及B1/secant lane。开发与确认合计48个case/lane/F状态，静止secant的3个状态因无观测速度缺失；其余45个实际输入均无接受planes，A/S在事件或plane阶段退出，没有真实搜索、无H运动可验、无可执行LOO。静止T仅3事件，保留事件不足。

开发选择[selection_development.json](real_frontend/selection_development.json)固定为none，保留F0基线，没有因F2训练簇较小残差强选。确认不重新选择。完整每seed latest、cluster、空间比、rank、RMS/gate与退出原因已保存；开发5046行支持重建见[逐seed CSV](real_frontend/development/seed_diagnostics.csv)。run1 F1中位latest11、F2簇7；run2分别6/4；平移分别3/3，剩下达到10点的簇也未通过原最终残差线。这不是仅统计有效plane数量。

有上限的27组合额外诊断固定原3个非静止开发窗及相同原T物理seed，比较T、前6ms T∪H、full10ms×3F，只做提取/rank/LS诊断，不search、不改变选择。27组合均无接受planes，原T结果与主结果一致。[诊断结果](real_frontend/split_diagnostic/summary.json)保留来源pixel数、邻居/簇/残差。F1/F2因果且seed固定在前6ms，所以T∪H和full10ms相同属于预期；不能称检验了后4ms的新seed。该消融只说明在这些固定seed上恢复像素/时间支持仍未得到通过门限的流，不能泛称所有数据不存在flow或证明模型错误。

### 一次完整108窗回归

遵循[最终范围协议](real_frontend/protocol_confirmation_regression.md)，主B1 lane、完整common_ids和同T/H/V，全部3F×A_literal/S_comp=648方法状态，另108个共享G训练参考。8个已算B1窗经radius/source SHA和输入hash核对后复用，不重算；secant只限固定8窗。回归只global fit/H/V，LOO明确not_evaluated，计数/资格为null，没有当0或通过。

|648状态类别|数量|
|---|---:|
|geometry_unavailable|120|
|T事件不足|84|
|plane不足/初始化不可用|444|
|A/S可用实录初始化或实录搜索|0|

F0/F1均0planes。F2只有5窗共7planes，最多2个、最高rank2：run1 4.365s 1个，run3 2.625s 2个、2.695s 1个，run4 2.425s 2个、4.495s 1个。其H总预测覆盖仅1.7–6.3%、区域支撑最多2，远未达预定覆盖/满秩要求；不能凭被预测H的conditional低残差称运动可信，亦未把回归中个别非空结果升级为winner。

G参考38窗native有效、39窗通过共同质量。差出的run5 4.025s被原G多初值竞争门限拒绝，仍保留拒绝而不宣传有效；这与A/S没有执行多初值检查的状态明确分开。本轮T最多3719事件，未触5000内部cap，三路线优化输入规模一致。该G输入为前6ms训练子集，不能与旧10ms结果逐项声称不变。

## 适用限制与最终验证

H/V事件隔离条件于旧冻结几何；secant还用窗末几何信息，因此不是完全因果在线预测或独立几何真值。H/V无可验运动不是测准零，不是误报改善；纯平移不能被translation-only稳定解释也使负例“成功”不可识别。A/S在进入未来线之前就因前端/秩失败，不能用V门限解释allzero。

[regression verification](real_frontend/regression108/verification.json)已通过：648/108完整；264个有输入case/F的T/H/V集合互斥且并回完整common_ids，所有拟合支持仅T；复用8窗radius、source SHA及输入内容一致；全部A/S失败无搜索/G回退；LOO为not_evaluated/null；参数/选择hash、task13封存及88个输入来源未变。此前[合成](real_frontend/synthetic/verification.json)、[开发](real_frontend/development/verification.json)、[确认](real_frontend/confirmation/verification.json)验证覆盖F0原接口一致、latest因果、H/V目标相同和几何/双向预测路径。

**结论：task1完整正例成功；task3语义边界明确；task2仍未取得实录可验证角速度。** 本轮不支持升级F1/F2、不支持Sony官方精度复现或实录提速声明。绝对轴/RPM精度、零速置信界和在线性能依然不可识别。三任务的代码、失败、测量诊断和完整紧凑结果均交付；全部既定计算结束，不再追加实验、scene、前端、图表或参数搜索。
