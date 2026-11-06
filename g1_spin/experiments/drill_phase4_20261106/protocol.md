# 阶段4 task1+3预定协议，2026-11-06

GPT-6 Astra xhigh审阅，GPT-6.1 Sol high执行。先仅完成合成几何事件全链正例与论文语义核查，再交Astra阶段审核决定实录前端；不提前运行实录变体。仅写本目录；阶段1/2/3、生产源码、原始归档不改。引用阶段3flow/methods冻结副本与hash，不改原工程邻域、平面门限、初始化可用条件、IWE/搜索/质量门限，不调参凑通过。现有venv正常执行，不加依赖、不委派、不问用户。

## task1：同一整数events贯通提取、LS与搜索

几何edge crossing生成器仅使用独立SciPy Rotation与根求解，不调用估计器。r120px，c=(640,360)，时间[0,10ms]，旋转参考5ms。材料q(t)=R[−omega(t−5ms)]lift((pixel−c(t))/r)，visible z<0。六个材料弧锚点xy=(.6,0),(0,.6),(−.5,.2),(.2,−.5),(.4,.4),(−.3,−.4)，image normals=(1,1),(1,−.3),(.2,1),(.7,−1),(1,.1),(.4,1)。用lift Jacobian J构造n3=normalize[J(JᵀJ)^−1 n2]，每edge仅保留距anchor≤.14rad材料弧。对整数像素解n3·q(t)=0，固定25us扫描及brentq，检查根/弧/可见性；输出整数x/y、四舍五入integer us和亮度跨越极性。固定pixel枚举以投影弧包围盒加2px，不根据输出扩大patch。标签edgeID/omega等只写evaluator truth，runner仅接events及明确几何；生成器不依估计结果选事件。此为理想几何边缘crossing，不是完整DVS传感器仿真。

固定场景：omega=(10,20,30)、反向、第二omega=(−20,30,10)；同texture零旋转（无event应无测量）；纯平移cdot=(400,−300)px/s；baseomega+同平移；旧169单法向精确平面作为秩退化控制。固定truth holdout：整数像素按((x+2y)%5==0)分出20%留出事件，训练只提取其余events，所有搜索也只用训练events，evaluator使用未参与提取/优化的完整holdout事件检验局部plane预测与全局共同omega的crossing约束。固定split不依输出改变；记录物理edge覆盖仅evaluation使用。对照完整事件extract另作生成正例支撑核验，不将truth边缘ID输给前端。

三旋转正例使用冻结extract→A_literal/S LS→有向axis scan+实际1次NM，G单独参考，无trueinit或G回退。检查S init axis≤10°且speed相对误差≤15%，最终≤5°且≤10%；反向有向符号正确。A_literal至少base/反向实际搜索，偏差如实记录，不要求等于S或truth。上述仅为合成实现验收门槛，不是实录精度。若固定生成器不通过，先核根/单位/约束/搜索，立即报告失败层，不无限改参数或降门槛。

## task3：语义、平移与约束

有限一次核查主文v2 §3.3.1、官方公开pseudo_code和可得补充，记录实际公开范围。阶段3 A原语义保留并明确名A_literal，不静默改变。保留原非等价反例，不宣称Sony未公开实现错误。

raw aperture约束r(s×[a,b,0])·omega=1−gᵀcdot。验证直接raw梯度修RHS与整个events坐标先减c(t)再extract两路径；不得用f_normal−cdot冒充完整2D flow。A_literal只运行在统一去平移并重新提取的events上，不自造A速度补偿。公式、前端和几何速度来源单列。

解析语义对照：h=(I−ssᵀ)r[gx,gy,0]，b=1−gᵀcdot，则hᵀsdot_rot=b，最小范数切向normal velocity=b*h/(hᵀh)，P=h hᵀ/(hᵀh)。逐条对应标量约束，但直接stack投影LS引入1/||h||残差权重，不能称整个LS等价；保留b=0，不除b构造法向。只做解析对照，不新增第4实录算法。

## 核验和停止

保存协议/生成配置/引用hash，原始events和truth分离，生成器根误差/整数时间误差/边缘覆盖/留出一致性、rank/condition、LS和最终误差、实际NM计数、搜索scan/zero及原候选。至少一例独立重建LS/球面xyz/目标值，失败不回退。合成平移路径与零/秩退化控制完整保留。源文件和旧结果hash不变。task1+3完成或发现固定正例无法满足时立即向父代理/Astra报告，不开始实录前端，不上传完整RAW或超大文件，不动源树无关main改动。
