# 开发none后的固定确认及108回归（Astra批准的最终范围）

开发selection_development.json已冻结selected_frontend=null，不升级F1/F2。27个同T物理seed的T/T∪H/full10ms×3F测量诊断完成后，不用于重选。四个确认窗仍按原三个前端及B1/secant两lane运行；若初始化可用则执行既有search/H/V和真正留区重拟合，失败保留，确认不调参。

108一次回归覆盖全部原预定window，包括geometry_missing：主B1 lane，完整common_ids，同预定T/H/V。三F×A_literal/S_comp合计648方法状态，另108个共享G训练参考。已完成8个开发/确认窗的B1结果，核对完整输入/几何/参数hash后复用，不重算，不称独立测试。secant只限固定8窗，不扩到108。

回归只global extraction/LS/search/H/V，LOO明确not_evaluated，其计数和资格标null，不能记0或通过。可用新candidate仍只作诊断，不重选winner。所有geometry_missing、T不足、plane/rank/质量失败和运动不可评都保留。G同T/lane仅参考，绝不作为A/S回退。

H/V都条件于阶段1/2/3旧冻结几何，几何可能使用整个估计窗信息；secant甚至使用窗末观测。因此这是给定几何假设下的事件T/H/V隔离，不是完全因果在线预测或独立geometry验收。旧108仍来自同五段数据，不能当新独立测试。已知真实合成平移未过V工程线说明指标适用边界，不能据此证伪真实运动；前端allzero发生在V检查之前。

全部规则仍以protocol.md顶部实录前审阅修订和parameter_freeze.json代码为准。禁止新增第4F、scene、阈值搜索或扩窗。完成确认与一次108后整合三任务交Astra终审，不自动开启后续实验。
