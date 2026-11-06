# 阶段5原始结构检查点（完成2026-11-07）

**旧F0的主要拒绝并不只是线性平面优化不足：多数邻域里，同一pixel、同一polarity在多个时间反复发出event，任何单值时间函数都已不可能达到原RMS门限。** 这一结论仅解释保留重复event的F0，不能解释去重复后F1仍失败，也不等于识别出真实运动、判定硬件错误或证明所有事件无flow。

本轮保留`drill_phase5_20261106`启动标识，实际执行完成日期为2026-11-07。Astra xhigh审阅方案，Sol high执行。只做原始观察、variance下界和测量正例校准，没有CMax、阈值放宽、新模型拟合或108重跑；确认窗未分析。

## 固定选格和原始事件

3个开发窗沿用stage4 B1/r/完整common_ids：run1 1.495s、run2 2.095s、run5 4.095s，静止活动控制run5 1.525s。只用前6ms event数，在绝对sensor12×12px格按count降序、格坐标固定排序，≥20events、中心距≥16px、最多6格。最终选6/6/3格，静止0格；空槽和不足保留，没有按残差/omega或后4ms选格。

选定格后，从完整官方解码按精确pixel/time取该格全部10ms原始event及物理ID，另标是否属于旧common支持。每格保存极性分离x/y/t、20个0.5ms切片、pixel计数/时间raster、first/median/latest时间面和同pixel间隔，见[选择表](patch_selection.json)、[全部量化](patch_metrics.json)、[紧凑CSV](patch_summary.csv)与`patches/`。后4ms图像仅作原始观察，不作后续训练/选择，新的holdout评价必须按另行冻结协议执行。

## 不可消去残差的直接证据

旧gate为RMS≤max(20us,0.15×std_t)。非20us floor分支要求R²≥0.9775。对相同pixel同pol分组，有total variance=within+between；同组所有event只能共享一个预测时间，最优单值函数也至少留下sqrt(within)的RMS。它比线性平面更宽松，因此该下界超过gate时，换更充分的线性优化也无法通过。

独立从stage4保存的原F0物理邻域ID重建结果：

|开发窗|≥10点邻域|单值函数下界已超原gate|下界/gate中位|
|---|---:|---:|---:|
|run1 1.495s|269|256|3.9114|
|run2 2.095s|186|173|3.3404|
|run5 4.095s|95|71|3.5290|

执行者值见[逐邻域CSV](old_F0_neighbor_lower_bounds.csv)及[分组JSON](old_F0_neighbor_lower_bound_summary.json)。Astra独立低层审阅得到相同计数及约3.91/3.34/3.53比值，来源另列[reviewer_source.json](reviewer_source.json)，未伪称同一个getter返回就是独立验证。

15个按事件数选出的格，在前6ms的pixel/pol分组下界也全部超原gate：run1各格ratio约2.80–4.64，run2 2.57–3.34，平移3.29–5.01。多event组占event比例分别约77.3–93.7%、43.4–66.2%、35.7–60.0%；不应把first/median/latest时间图当作完整原始event的无损时间表示。

## 可见局部结构与不能确定的原因

- 同pixel同pol重复主要有毫秒级分离：15格的前6ms同组间隔中位约0.742–3.225ms，≤50us的间隔比例均为0。因而“主要由≤50us同pixel密集短burst导致”不被这些格支持；仍不能排除较长的多阈值响应、不同边缘/纹理多次经过、运动/曲率和照明耦合。
- 原始x/y/t和raster中，同空间格存在跨时间层的events，既有局部可辨走向也有与其不共面的点。例如run1 patch0两极性pixel重合Jaccard约0.58，而patch1/4/5的异极性空间最近距离中位约5.83/5.05/3.08px；run2前4格约1px且重合约0.19–0.22，平移3格前6ms均只有一种极性。它们说明不能把所有原邻域默认当同一条薄边缘同一次通过，但尚未进行edge归属或局部速度拟合。
- 同时空支持上的多值t(x,y)下界不能由曲面拟合消除；但这不区分“一条边缘连续对比触发多event”与“多条边缘/纹理轨迹”。去重复后的F1可能仍受时间选择、边缘混合、曲率/变速、空间rank和噪声影响，本检查点尚未证明哪一个是主因。
- 球外对照取距B1>1.5r的原始event。球内/球外50us bin率相关约−0.067、+0.029、+0.110，静止约+0.014；各patch与球外相关绝对值最大约0.33，没有观测到统一的强同步活动模式。球外含钻杆、手和背景运动，相关性或无相关性都不是独立照明/同步真值，不能据此排除局部光照影响。[完整活动序列](outside_activity.json)保留。

代表原始图：[run1 patch0 x/y/t](patches/run1_01495000_patch0_xyt.png)、[同patch raster/间隔](patches/run1_01495000_patch0_raster_intervals.png)、[run2 patch0](patches/run2_02095000_patch0_xyt.png)、[平移patch0](patches/run5_04095000_patch0_xyt.png)。所有15格和空极性均保留，没有为报告挑漂亮的格或丢弃失败。

## 测量正例与独立验证

[measurement_calibration.json](measurement_calibration.json)仅校准测量代数，不引入新的real门限：169像素精确平面回收原梯度，plane RMS约5.5e−13us、下界0；同pixel人为三个±200us重复的507events下界精确为163.299us，超过其旧gate104.103us；分开polarity下界0而错误合并为500us；20us floor分支正确保留。没有以这些控制设计新的空间模型资格。

[verification.json](verification.json)通过：15个已选格从完整NPY用独立整段布尔筛选复核全部物理ID/像素/时间/极性；独立计数排序复核选择只用前6ms；独立逐组循环重建within/between分解，和主实现一致；旧输入及来源hash未变。前次输出只出现NumPy bool不能序列化的落盘问题，已显式转换为Python bool后恢复同一固定数据统计，保留已绘图，不改变选择/数学/输入。

本检查点结论是：**完整重复events与严格单值时间面gate的结构冲突已得到直接证据；F1去重复仍失败的原因和可验证局部运动模型尚未解决。** 到此交Astra原始诊断检查点。后续R_all/R_latest×M_time/M_spatial的2×2表示/拟合及新holdout校准门槛仍待完整协议，未提前执行；确认窗、CMax和108均未运行。代码、统计、小型原始裁剪和图保留，完整RAW不上传。
