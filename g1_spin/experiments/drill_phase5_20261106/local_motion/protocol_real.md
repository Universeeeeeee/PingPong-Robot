# 真实2×2开发协议与资格（2026-11-07）

本开发仅固定3运动窗+静止活动控制，6/6/3/0原始已选格；每格3段×两polarity，360条四组合记录，失败不补格。确认run3/run4/run5三窗本检查点未分析；不做108或CMax。tau由已审clean oracle固定0.5px，污染/静态控制不抬tau。

T先按前6ms及80%固定rawpixel分出，H为其他20%pixel；后4ms为V。trainer仅接独立T NPZ，32次整rawpixel簇有放回bootstrap在representation选好之后执行，两fit使用相同抽样。每draw只需代数可解，不要求再有12 unique；≥24有效、原数据≥12T unique、时间可辨、normal方向P95≤15°、finite，才是T合格。n对齐原n时同步翻vn。delta_v为bootstrap |v_boot−v_hat|的P95，区间保存。

T motion标签：|v_hat|>max(50,2delta_v)为candidate_nonzero；|v_hat|≤50且区间含0为static-compatible；其他unidentifiable。S约束保留实际v_hat，不将static标签强制赋零。最小范数normal velocity为vn*n；vn=0时A_literal投影缺失，不能除vn反推g。

同T表示支持独立拟合最佳固定空间line作为static对照。H当前段至少4不同pixel，下一1ms future至少4不同pixel且2个半msbin；所有(pixel,pol,halfmsbin)内全event距离平均、单位等权。moving两边P95≤tau，且gain=(L_static−L_moving)/max(L_static,tau)分别≥.20，L为等权mean绝对法向距离。预测1ms位移|vn|*.001<tau时保守标分辨率不足。static低残差只能记stationary_constraint。

19次固定时间置乱在当前Tseg/H/短future各自内部进行，保留xy/pol及event集合；R_latest在置乱T重选latest，重拟moving/static。statistic=min(Hgain,Fgain)，至少15可评且真实statistic严格大于全部可评置乱值，才过permutation。它仍不能绕过T nonzero标签、距离/支撑和位移分辨率。

全局系统只用**全部T合格**约束，H/future/permutation绝不筛约束或重新拟共同omega。每组合至少12合格obs、4独立spatialpatch、rank3/cond≤100，leavepatchout仍rank3才允许试共同运动；同组nᵀu=vn拟2D translation并保rank/cond/LOPO。预测同heldpatch T末1ms模板及预定完整common支持V目标，固定patch bbox+15px，8个半ms片段、双向3px截断距离和pixel/pol/bin等权，比较omega/translation/zero。该额外预测只诊断，不因相对好一点覆盖局部验证失败或宣称新metric已经校准。

开发组合必须两个旋转窗均T充分、至少3patch有H/future/permutation支持的非零局部运动并有共同LOPO证据；负例支撑不足不是成功。若多个合格，比两窗中较弱gain，1e−6并列按Rall/time、Rall/spatial、Rlatest/time、Rlatest/spatial；无合格则none，不能改参数求通过。确认不重选。

时序证据：pre_real_freeze→development_training_freeze→development_qualification_freeze先保存T模型/bootstrap/全部资格约束及全局系统hash，再real_evaluation_freeze，最后打开H/future。原开发T源码冻结保留。Astra审阅发现未来sphere-invalid分支可能出现矩阵/LOPOmask行数不一致，后续源码已改为整系统geometry_invalid拒绝，保存旧源码及geometry_path_amendment；本dev无非法point，数值结果未重跑，旧freeze未覆盖。
