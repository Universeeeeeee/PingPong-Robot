# A语义、法向协向量与平移：保留原解释，不静默改算法

本次有限来源核查日期2026-11-06。主文v2 §3.3.1给出plane-fit二维flow、按r归一化、球面切向lift、P投影、LS、幅值搜索及NM流程，但没有写明时间平面梯度g如何成为所用三维normal-flow向量或其度量。因此阶段3的解释继续明确命名 **A_literal**：f=g/(gᵀg)，v=Jf/r，P=vvᵀ/(vᵀv)。阶段3原源与本阶段冻结副本均未修改。

来源范围：

|来源|实际可确定内容|不能据此确定|
|---|---|---|
|[主文arXiv v2 3.3.1](https://arxiv.org/html/2606.26780v2#S3.SS3)|二维flow→半径归一化→切向lift→P→LS→幅值搜索→轴/幅值细化|g→v细式、normal-flow度量、具体邻域/聚类/时间表面筛选|
|[官方ace_public approximate pseudo code](https://sonyresearch.github.io/ace_public/pseudo_code/)|每极性最近时间戳时间表面、CNN初值、球面CMax；页面声明仅近似关键流程|offline OF/P细式；不能把该时间表面说成官方offline front-end规定|
|[Nature 23页Supplementary Information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-026-10338-5/MediaObjects/41586_2026_10338_MOESM1_ESM.pdf)|RL、FAOC、serve及系统评估等；全文搜索optical/contrast没有匹配|未提供可解normal-flow歧义的offline flow推导|
|Astra对19页v2 PDF的独立只读审阅（父代理转述）|PDF不包含文中引用的Appendix A补充|未获得的补充内容仍不能臆测|

没有复制第三方实现。本次只检查上述来源，不泛搜、不反复追索。资料不足以确认Sony未公开offline实现究竟选何解释，更不能根据A_literal与S的差异宣称其实现错误。

对单位球s=(x,y,z<0)，image plane为t=a*u+b*v+c，g=(a,b)是时间对像素的协向量，单位s/px。完整图像速度为r*(omega×s)的xy分量加球心速度cdot，因此：

`r*(s×[gx,gy,0])·omega = 1−gᵀcdot = b_rhs`。

这里b_rhs是aperture约束右端，与时间平面的y斜率符号无关。raw事件提取梯度时保留该右端，不能随意令1。每个raw平面代表点使用其空间质心和拟合质心时间对应的几何中心位置。

把r*[gx,gy,0]投影到球面切空间，得到h=(I−ssᵀ)r*[gx,gy,0]。因为sdot_rot=omega×s且sᵀsdot_rot=0，故hᵀsdot_rot=b_rhs；该约束的最小欧氏范数切向速度为v_dual=b_rhs*h/(hᵀh)，方向投影为P_dual=h hᵀ/(hᵀh)。P_dual用h定义，**即使b_rhs=0仍有效**，不得先除b_rhs构造法向。

逐条约束与S相同，但stack `P_dual*(-[s]cross)`和`v_dual`时，其残差平方为`(S_row·omega−b_rhs)^2/||h||²`，即按1/||h||缩放标量残差。因而其LS等价于显式加权S，通常不等于未加权S。这里只作解析语义对照，未作为第4实录算法，也未替换A_literal。

A_literal使用Jg/(r*gᵀg)的方向；正确切向协向量的度量还涉及JᵀJ，二者一般不同。阶段3单点反例和六点LS反例保持原值，测试文件[semantic_verification.json](semantic_verification.json)明确核对原单点差异范数0.245973675。不能将A_literal改成dual后仍称原A。

平移处理明确分成两条路径：raw梯度修约束右端，或将全部events坐标减去cdot*(t−t_ref)后重新提取平面。对理想局部平面、b_rhs非零，变换后梯度是g/b_rhs；有限曲率、采样和支持变化使两条实测LS结果只应近似一致。b_rhs=0时去平移坐标上的运动边缘可没有可辨识t(x,y)平面，不强行除零。纯平移控制下，这并不等于测得准确零spin。

`f_normal−cdot`不能作为完整运动去平移：即使纯平移gᵀcdot=1，f_normal与cdot也一般只在同一法向投影上一致，其差还含非零切向分量。解析测试保留了该非零反例。A_literal在有平移的scene中**只运行统一坐标补偿后重新提取的events**，没有人为对其normal-flow速度做减法。raw RHS S与统一补偿S分别报告，几何速度输入为合成给定轨迹，不冒充实录检测真值。
