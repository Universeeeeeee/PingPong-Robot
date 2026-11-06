# Task2开发检查点：未找到合格前端，确认和108回归尚未运行

**F0全T邻域、F1因果latest time surface、F2稳健surface簇，在两个预定开发旋转窗均未产生通过原门限的local plane；没有前端可被选为合格方案。** secant诊断lane结果相同。A_literal/S_comp没有可用实录初始化或真实搜索，H运动和跨区域共同omega均不可评；这说明本次固定输入/前端/筛选组合没有提供可验运动，不能据此证明运动模型本身错误。

选择冻结为`selected_frontend=null`，保留F0作为既有基线，F1/F2仅为未通过的诊断候选。没有按planes数量或残差最小强选。四个确认窗仅准备了物理输入和几何hash，**尚未执行任何确认估计**；108的范围等待Astra检查点决定，没有擅自缩成旧行为或提前运行。

## 实录之前的测量修订和合成检查

[protocol.md](protocol.md)顶部“实录前审阅修订”是实际口径，覆盖初版正文相冲突的汇总描述。初版仅六scene合成输出留在`synthetic_prereview_v1/`，代码/协议/freeze快照留在`pre_real_review_v1/`。Astra审阅后，在首次实录运行前修正：latest必须先在全部T因果history按sensor pixel取值、再查lane半径；source出球体明确geometry_invalid；背面预测保留3px惩罚；双向距离独立资格与等权mean Chamfer；spin同时优于zero和translation-only；被拒绝LOO不算通过；H增加固定总覆盖/三区支撑要求。[六项实际路径测试](path_verification.json)通过，新[parameter_freeze.json](parameter_freeze.json)随后冻结，未以任何实录表现改变参数。

保持原task1六场景，使用本次严格T/H/V（不是task13空间20%split）。36个scene/lane/frontend状态已完成。三个旋转正例×三前端的9个S案例全部贯通LS+search、通过原合成轴/幅值门槛和双方向V线。base/reverse/second的F1、F2分别接受197/185/233个planes，H法向距离中位约0.0107/0.0104/0.0107px，H总预测覆盖83.8%/83.3%/94.0%；F0为277/272/291planes，仍与原冻结extract在同T上逐项一致。它们证明新前端在干净事件边缘上可执行，不表示实录有效。

纯平移合成控制需单列边界：secant translation-only未来双方向中位1.793/1.587px，forward相对zero改善19.81%，未过各向≤1.5px/≥20%线；固定1ms源模板和离散稀疏crossing在这个慢平移场景未被该工程线验收。没有因此降门限、改模板或宣称补偿后0plane就是平移成功。

## 完整物理支持与T/H/V

从stage3 `common_ids`取完整事件，保留原物理ID；不使用旧5000截断结果。先split再surface/cluster：T为前6ms且不在固定2px棋盘H组；H为同前6ms其他pixel；V为末4ms全部。seed/邻居/triplet/LS/CMax只读T。诊断secant速度事先从旧几何窗首/尾求得并冻结，不能当真实中心速度。

|开发窗|完整events|T / H / V|secant速度px/s|主lane各F接受planes|
|---|---:|---:|---:|---:|
|run1 1.495s|5460|1709 / 1714 / 2037|(−44.64,−73.12)|0 / 0 / 0|
|run2 2.095s|2809|849 / 824 / 1136|(−59.17,−56.35)|0 / 0 / 0|
|run5 4.095s平移|762|241 / 229 / 292|(−258.06,−37.24)|0 / 0 / 0|
|run5 1.525s静止先验|18|3 / 5 / 10|无观测secant，留缺失|0 / 0 / 0；事件不足|

24个预定case/lane/frontend状态中21个有lane输入，静止secant的3个前端状态均明确缺失。A/S共42个实际route状态均在事件/plane阶段退出；没有失败回退G、没有真实A/S NM、没有伪造H/LOO预测。raw S_RHS语义诊断同样0接受planes；主A_literal/S_comp共用相同坐标planes，诊断lane先统一去平移再重新提取。

## 失败发生在哪一层

保存的每个seed记录包含raw邻居、因果latest邻居、cluster保留/丢弃物理ID、最多64triplet的可用数、空间支撑、拟合rank、时间std、残差/门限。独立从保存的T支持重建了[5046行逐seed诊断](development/seed_diagnostics.csv)，没有重跑新窗口或重新选择cluster。分组见[seed_diagnostics_summary.json](development/seed_diagnostics_summary.json)。

|主B1窗 / 前端|raw邻居中位|latest邻居中位|最终簇中位|≥10邻居/簇的seed|仍因最终残差退出|
|---|---:|---:|---:|---:|---:|
|run1 F0|28|28|28|269/300|269|
|run1 F1|28|11|11|176/300|176|
|run1 F2|28|11|7|73/300|73|
|run2 F0|13|13|13|186/300|186|
|run2 F1|13|6|6|46/300|46|
|run2 F2|13|6|4|1/300|1|
|平移 F0|7|7|7|95/241|95|
|平移 F1|7|3|3|7/241|7|
|平移 F2|7|3|3|0/241|0|

F1因果latest和F2聚类进一步减少了支持。剩下达到≥10的簇仍未通过原RMS门限，F2在run2只有一个足够大的簇，也不通过。secant lane完整计数均保留，未按表现择lane。T/H/V分离确实减少训练事件和邻域，但本次没有增加未分离F1/F2对照去唯一归因；阶段3同窗完整旧输入的F0本来也是0plane，不能简单把失败全部归因split。更不能以训练簇残差有所降低就称获得了运动。

## 留出与负例

因没有通过前端筛选的planes，所有H distance为不可评、总H预测覆盖0，没有把它记成距离0；全局LS/search及真正留区重提取都不可执行，所有LOO成功计数为0而非“验证通过”。这正是资格失败原因，不允许绕过前端给旋转候选补证据。

zero和translation-only仍按同模板、八个V片段、同ROI和同目标ID评价。run5开发平移的secant translation-only双方向中位为2.950/2.479px，forward比zero仅改善1.67%，reverse反而约恶化10.88%；没有通过“平移可解释V”的工程线，故不能把前端无flow当负例成功。静止源模板不足，未来不可评。两个旋转窗的secant-only预测也未过线，这些速度只是旧几何假设。

同T/同lane的G仅作单列参考：run1两个lane有被接受候选；run2两个lane在原质量门限下拒绝；平移拒绝、静止事件不足。T前6ms/完整common支持不同于旧10ms/最多5000输入，这些G状态不能冒称旧108完全一致回归，也没有用于选择F1/F2或作为A/S回退。

未来资格实际要求**两个方向各自**median≤1.5px且各自优于zero≥20%，并同时优于当前lane translation-only；报告等权双向mean Chamfer及各方向覆盖，拼接median仅诊断。source非法球面不能clamp，背面source不能删除。前端选择尚未涉及V筛“好目标”，每个模型和两个baseline目标物理IDs已独立核查相同。

## 选择与核验

[selection_development.json](selection_development.json)将没有合格前端的事实及来源hash冻结；两个开发旋转窗都未达LS/search条件，因此所有后续资格均未通过。F0保留为原基线，不把F2更小训练簇残差当选胜理由。确认和108下一步应由Astra明确范围：固定确认仍应保留三个前端的表征/邻域证据，108用作同录制历史输入回归而非新独立测试；本检查点未提前执行或缩减。

[synthetic verification](synthetic/verification.json)及[development verification](development/verification.json)通过：F0同T结果与原extract一致；全部seed/surface/cluster只含T，latest全history先于半径；H/V无拟合泄漏；预测模型/zero/trans同target物理IDs；拒绝measurement字段null；已封task13的61个文件及旧输入来源hash未变。所有输出都是局部边缘运动假设诊断，0plane不证明G1不支持、Sony失败、真实无运动或模型本身错误。

完整压缩detail、训练簇与丢弃IDs、plane数组、候选和计时保存在`synthetic/`和`development/`。实录开发只到既定4窗，参数未调、未扩窗、确认未运行；到此交Astra检查点。
