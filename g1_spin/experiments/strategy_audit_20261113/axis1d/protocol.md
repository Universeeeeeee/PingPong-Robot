# 已知轴一维 CMAX

使用第二次电钻记录的固定轴/反转/姿态窗口，使用原阶段3冻结G输入；不改事件、ROI、中心、半径、最大RPM或原生产估计器。拟合只搜索一个有符号标量s，omega=s*a。保留搜索中的全部局部峰。

axis来源分开报告：run1_vertical_image/run4_vertical_image是严格粗物理假设[0,±1,0]，因为录制标签称钻轴画面竖直，但未给出深度分量；run1_free_mean/run4_free_reverse、run2_data_mean/run3_data_mean由此前自由CMAX有效候选汇总，只能诊断一维幅值稳定性，不是独立轴真值。

与自由三维CMax逐窗比较RPM、符号、局部峰和相对zero对比度；不以247RPM弱参考选峰、不以自由结果选择粗轴。该实验检验可行性和低维条件性，不验证未知轴的三维精度。失败保留，历史结果不覆盖。
