# Task2有限前端协议：合成检查→开发选择冻结→检查点→确认/108回归

GPT-6 Astra xhigh审阅方案、GPT-6.1 Sol high执行。当前只授权合成控制检查及四个开发窗；开发选择文件写出后交Astra，再按下发要求执行确认与108回归。仅写本子目录，不改task1+3、阶段1/2/3、生产或归档，不新增依赖、算法路线或窗口，不用实录反复调阈值。

## 固定输入与split

开发：run1 1.495s、run2 2.095s、run5 4.095s、run5静止1.525s。确认仅预备输入hash，不求结果：run1 1.425s、run3 1.495s、run4 2.495s、run5 5.095s。使用stage3固定B1/r和完整common_ids物理支持，不是最多5000截断。时间均相对该录制首CD。保持原物理event IDs、pixel、polarity和时间来源。

先split再提取：T=[0,6ms)且不在hold像素组；H=[0,6ms)的其余；V=[6,10ms)全部。固定2px棋盘组为`(floor(raw_x/2)+floor(raw_y/2))%2==0`分给H；T用相反组。F0/F1/F2的seed、邻居、surface、triplet、簇拟合、LS和CMax均只读T。所有T事件用于提取和搜索，不额外下采样；既有cmax若T>5000会原样内部确定性截断，记录其实际优化IDs，不用H/V补入。该T/H/V是未来预测评价，与task1空间20%holdout不同，也不作为独立录制测试。

## 固定三前端

F0：冻结全事件邻域4px、±2ms、同极性、最多300均匀T seed。

F1：同极性因果surface，先限定≤seed时间且age≤2ms，每个原物理pixel保留latest T event；同一时刻并列按较大物理ID确定，半径4px在当前坐标lane查询。pixel identity永远用原整数pixel，不把去平移小数坐标当新pixel。

F2：仅在F1候选邻域上seed锚定稳健平面簇。枚举seed+两邻居的组合，稳定次序中均匀linspace选最多64triplet；xy三角面积≥2px²。候选法向像素距离`abs(dt−g·dxy)/||g||≤.75px`定义训练簇，以最大支持/最小median全邻域距离/固定候选序择一次，然后仅一次LS重拟合。该cluster不在H/V上挑inlier。保存训练邻域、保留簇、被弃T事件和退出原因。

三者最终均保持原≥10neighbors、二维支持比≥.08、平面秩3、时间std≥1us、RMS≤max(20us,.15std)、梯度≥1e−8及0.9r范围等筛选。共同LS原≥12planes/rank3/condition≤100/有限非零初轴，保留原S native残差；公式/质量门限不变。没有第4frontend或失败回退G。

## 两几何lane与公式账本

主lane B1：固定stage3中心，vc=0，A_literal/S_comp共用raw坐标planes。

诊断lane secant：vc=(既有几何end−start)/10ms，先冻结来源和值；中心假设c(t)=B1+vc*(t−5ms)，只为假设，不按spin分数选择。缺少实际插值支持则诊断lane缺失，静止人工位置先验不伪造速度。A_literal和S_comp共用统一去平移后重提取planes，同lane G也用同坐标及T IDs。为同时保留同坐标公式对照和raw平移语义，另存原始planes上S_raw_RHS，右端1−gᵀvc，球面点用c(t_local)。它是同lane独立语义诊断，不参与前端择胜，不是第4frontend；这明确区别于原计划将rawRHS单作主诊断路线。A_literal绝不直接f_normal−vc。

## H局部运动验证

对每个T提取plane，用固定4px距离、同极性、相对seed时间±2ms先选该plane附近全部H（不因plane残差删H），评价法向像素距离`abs(time residual)/||g||`，记录分布、中位/P95、≤1.5px覆盖及H的全部物理IDs。同一H可被多个plane预测，保存配对；另报按H event聚合的中位距离，避免只数训练簇或最好一次。没有H支持留不可评。

## V真正未来预测

仅T末1ms=[5,6ms)作为源边缘模板，每个源event按自身时间lift到球面，用完整`R(omega*(target_time−event_time))`传播到V八个0.5ms片段中心（6.25至9.75ms），再叠该lane中心轨迹。V既不进入IWE、拟合、选择peak、中心或对应，也不只挑预测好的target。

全局ROI在模型前固定为所有T raw像素bbox+15px、裁传感器范围；区域留一评价ROI为该区T源event raw bbox+15px（无源则不可评）。所有模型、zero和translation-only共用模板、target IDs与ROI。zero=(omega0,vc0)，translation-only=(omega0,当前lane vc)。按极性双向最近邻Chamfer，距离截断3px，给所有预测点和所有ROI内V目标评分，未匹配极性距离3px。每片至少10源event、10target event才可评；8片均列出，不从缺失/失败片填值。记录双向距离中位/均值/≤1.5px覆盖，整体用全部有效片距离汇总，中位≤1.5px且相对同ROI zero中位改善≥20%为工程预测线；至少4/8片可评才可对该线判通过，少于4片记不可评。平移与zero同模板，不能把补偿后无plane自动算负例成功。

## 空间区域和留区重拟合

在T raw像素相对B1中心的固定24px方格中计数，以T事件数降序、格坐标次序择最多6个区域；每区至少20个T event，贪心区域中心距离≥24px，计数选择只用T，不看H/V或omega。区域bbox为该24px格内实际T raw bbox。初始运动候选需要至少3个拥有接受plane seed的已选区域，且T plane normal的无向角至少两簇相差20°；不把同一区众多plane当独立证据。

每次留一区删除该区全部T events及bbox外扩4px缓冲，再从其余T重新提取/LS/search；被留区域的T末1ms源event模板只交evaluator，绝不用于重拟合。仅对全局有可用初始化/搜索的route执行其最多6次留区重拟合；其余直接报告不可评。每次rank不足/不足planes保留不可评，不缩缓冲/改区域追求rank。评价被留组固定ROI的V，记录逐区H/V及共同omega差异；组只是空间/法向分离候选，不宣称独立物理边缘或独立试验。

## 资格和开发选择

正例先在既有task1六scene上用同T/H/V检查F1/F2与F0，报告通过哪一层；若既有清洁正例都失败，应定位而非放门限。参数一次冻结后才做开发实录。合成成功仅实现前提，不等于实录通过。开发选前端只看主B1：两个开发旋转窗均有LS+search、≥3region/normal支撑，H中位≤1.5px且覆盖≥50%，V工程线通过，至少3个留区结果可评且通过V工程线；开发平移没有被未来支持的错误非零（允许无plane但此时不能当平移解释成功），静止不足保留。无candidate满足则明确不强选，F0保留基线、仅列诊断结果，提交Astra决定如何做一次固定108回归。

若多个合格：先比较两个开发旋转窗V相对zero改善的较小值（越大越好），并列按F0/F1/F2顺序；不按plane最多选择。选择一次写selection_development.json并冻结hash后停止检查点。确认旋转需≥2/3且必须含反转，两个平移不得出现被未来验证支持的错误非零；确认不再调参数/重选。所有失败和候选/measurement/native/common状态、包含flow的成本保留。108仅旧输入回归，不能当独立测试。
