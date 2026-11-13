# 前段拟合/后段预测与中心速度剖面

结论见上级report.md。主输入：9个先前已观察的实录窗＋2个理想点轨迹正例，前60%拟合omega、后40%评估。保留整窗冻结几何/ROI条件，非因果端到端验收。

结果：7个native有效实录旋转候选均改善原像素未来预测27.1%–52.6%；但原协议要求IWE余弦与像素指标双胜，正式通过4/9。drill2原门限拒绝，不能把其预测改善当有效测量；drill5纯平移拒绝候选只微小改善0.46%，不构成物理证据。27个真实训练置乱均no_rotation_evidence。三置乱不提供统计显著性。

两个真实中心速度剖面共10个假设，保留全部事件与数据划分。圆盘205.38–240.93RPM对应距离1.2335–1.2549px，电钻248.16–280.02RPM对应0.8149–0.8292px。人为速度偏差不是独立测得中心误差，不能把RPM范围当置信区间，也不按最低H误差选择模型。

核验：verify.py从保存输入重新重建球面、训练/留出，独立Rodrigues+brute-force cdist复算全部像素评分，与cKDTree实现一致。verify_profile.py同样检查剖面，所有事件保留。原reversed_training仅反点序，是非均匀时间下不成立的严格时间反转诊断，不用于结论；正确的t'=t_first+t_last-t_reverse、p'=p_reverse变换已独立验证warp等变。

协议写points/polarity打乱，但CMAX没有极性输入，实际只打乱point与time对应；评价始终使用真实训练模板和真实极性。这个口径差不改变CMAX输入，但必须明确保留。

代码入口run.py/verify.py/profile.py/verify_profile.py，运行估计的输出拒绝覆盖。results.json完整保留优化历史，summary.json为逐case裁决，profile_results.json为参数族，两个verification.json保存独立核验。预测几何真值仅synthetic；实录标签不升级。
