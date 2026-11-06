# 阶段5首检查点：原始局部结构与单值时间面下界

启动标识保留drill_phase5_20261106；实际执行/完成日期2026-11-07（Asia/Shanghai）。Astra xhigh审阅方案，Sol high执行。仅原始/数学测量，不运行CMax、不新增前端或阈值、不做108回归，不改旧证据、生产代码或RAW。

开发固定run1_01495000、run2_02095000、run5_04095000；确认run3_01495000、run4_02495000、run5_05095000预留不读取结果；run5_01525000静止活动控制。沿用stage4冻结B1/r和完整common_ids物理支持，不用旧5000截断。开发选择patch只用前6ms common支持事件数，不使用旧plane残差、omega或未来event。固定绝对sensor12×12px格，以count降序、格x/y顺序；每格前6ms≥20events，贪心中心距≥16px，最多6格。格不重叠，不足/空留记录，不扩格追结果。

每个已选格完整10ms原始event从官方整段归档精确像素/时间筛选保留物理ID，另标是否属于common_ids；这可避免对同pixel重复event再次mask。前6ms用于本检查点的variance等结构量化，6–10ms只可视诊断，绝不用于后续训练或选择；后续训练/holdout必须重新遵守单独冻结协议。

每patch保留按极性x/y/t散点、20个0.5ms原始event切片、按pixel/polarity次数和时间raster、first/median/latest时间面、同pixel同pol间隔。附原始event率50us bin、每组相邻事件间隔≤50us的时间连通簇（仅描述，不筛inlier）、正负极性pixel重合和最近邻空间距离。时间簇50us、时间统计/图参数固定，不按结果改。

数学诊断：旧gate RMS≤max(20us,.15×std_t)，非floor分支意味着R²≥1−.15²=.9775。对相同pixel同polarity分组，total variance=within-group+between-group；任何单值t(x,y,pol)最小可能MSE为within-group variance，RMS下界sqrt(within)。这比线性平面更宽松，若下界已大于原gate，则在该完整重复events上任何单值txy都不可能过门限；不是优化器不充分。分极性和联合pixel/pol两种账本都保存。它只能解释使用重复events的F0，不能解释已去重复F1仍失败。

独立重建stage4三个dev B1 F0的全部保存邻域（同T支持、同pol），复核≥10points下界/gate。审阅者转述的Astra值（run1 256/269、run2 173/186、平移71/95，中位ratio3.91/3.34/3.53）另标来源，与执行者独立计算分开，不能冒称同一层getter为独立验证。

球外对照取同完整10ms、距B1中心>1.5r的event，仅比较50us bin同步活动、极性及总率；含钻杆/手/背景运动，不作为独立照明/同步真值。不能从相关性单独归因为照明。静止窗仅活动对照，无≥20events格则保留0选格，不换时间。

停止于该固定原始结构/测量下界检查点，报告能支持/不能排除的burst、多edge、曲率变速解释及观测边界，不机械放宽gate。Astra后续2×2 representation/fit及新holdout校准方案尚未执行，不能提前取代当前诊断。输出小型裁剪/统计/图/源码，完整RAW不上传；原工作树main无关改动不触碰，后续按父代理安排更新PR4。
