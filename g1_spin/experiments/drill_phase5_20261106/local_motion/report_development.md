# 2×2开发检查点：训练方向有改善，尚无通过留出检验的局部运动

**已固定完成360条开发2×2记录，selection=none。** 直接pixel法向模型比逆时间模型更容易得到稳定训练方向，但所有组合仍无H/下一1ms future/置乱共同支持的非零运动，也无被验证的stationary constraint；不能把代数可辨或训练残差降低当作真实角速度。确认三窗暂未打开，CMax和108未运行。

合成oracle tau=0.5px gate先通过，随后把15个既定格分成独立T文件与H/V评价文件。训练脚本只读T，32pixelcluster bootstrap和全T资格/全局系统分别封hash，之后才评价H/future和19次各split内部时间置乱。支撑、stationary/nonzero、static对照和组合资格细则见[protocol_real.md](protocol_real.md)。H/future从未删训练约束或重拟共同omega。

|开发窗/组合|T-qualified局部数|有方向支撑的patch|T candidate-nonzero|最终非零motion支持|
|---|---:|---:|---:|---:|
|run1 Rall/Mtime|8|4|4|0|
|run1 Rall/Mspatial|14|5|6|0|
|run1 Rlatest/Mtime|5|4|4|0|
|run1 Rlatest/Mspatial|13|5|6|0|
|run2 Rall/Mtime|0|0|0|0|
|run2 Rall/Mspatial|8|4|0|0|
|run2 Rlatest/Mtime|1|1|0|0|
|run2 Rlatest/Mspatial|6|4|0|0|
|平移 Rall/Mtime、Mspatial|0 / 4|0 / 2|0 / 0|0 / 0|
|平移 Rlatest/Mtime、Mspatial|0 / 3|0 / 2|0 / 0|0 / 0|

每个运动窗分别36/36/18个位置-时间-极性观测，再乘四组合；静止原始18events没有≥20events格，保留活动不足，不能当零速成功。双极性与同录制分段都不是独立重复试验。

## 失败分层，不只报allzero

[development_failure_attribution.json](real/development_failure_attribution.json)对16个case/组合提供全部边际通过数与按固定顺序的首个失败数：model代数→T资格→T nonzero→H支撑→future支撑→H P95→future P95→H gain→future gain→permutation→位移分辨率。这样可区分支撑不足、距离过线、相对static没有改善、时间置乱可解释或位移太小。

T上run1有少数candidate_nonzero，但当前2ms单平面在H/短future上的P95明显超过0.5px；run2即使空间normal稳定，也没有T上可辨非零speed。平移只有少数稳定方向而速度不可辨，多个观测H/future空间支持不足。完整H/future P95分位、mean gain分位与19置乱统计均保留在[development_evaluation.json](real/development_evaluation.json)，不删除坏事件、不挑最佳时间对应。所有stationary标签也必须独立支撑，不能因小vn默认成功。

这比阶段4的“原time RMS全拒”更具体：表示与pixel拟合改善了部分训练normal，但在当前12px/2ms/1ms预测尺度下，**一个恒法向/恒法向速度局部层尚不能解释全部未拟合原始event**。多层触发、纹理/边缘混杂、曲率/变速或采样影响仍可能并存，不能据此宣布硬件无flow或证明某一种物理原因。

## 全T共同系统与公平heldpatch诊断

仅run1两个M_spatial组合满足≥12 T-qualified、≥4patch、rank3/cond≤100及LOPO仍rank3。Rall为14obs，ω诊断约(16.69,−14.29,−17.84)rad/s、T aperture RMS498.3px/s；Rlatest为13obs，ω约(−4.95,−32.05,−2.75)、RMS506.6px/s。数值只是全部T约束LS结果，不是被接受物理measurement。

同一T约束还拟合2D translation，每个heldpatch用其自身T末1ms作仅evaluator模板，其他patch约束重拟的omega与translation对**同一固定ROI和V target IDs**进行完整传播/平移、zero对照及3px截断双向距离；10个heldpatch结果见[development_global_LOPO_predictions.json](real/development_global_LOPO_predictions.json)。没有按模型选bbox、目标或对应，也没有用heldpatch参与其全局拟合。各解释有时相对距离不同，但局部留出非零证据为0，因此不为omega作额外rotation解释，不强行进入A_literal/CMax链。

run2所有组合均不足12 T-qualified，已经不满足开发两个旋转窗共同资格。故[selection_development.json](real/selection_development.json)严格冻结none；没有借H筛出几个好planes再拟omega，也没有按训练normal最多选择Mspatial。负例支撑不足不计成功。

## 独立核查与代码修订

[verification_development.json](real/verification_development.json)通过：360模型/16系统完整、两fit同支持、T/H/V物理ID互斥，912次等权mean/P95独立重算一致；32bootstrap和19permutation数目完整；全局约束集合严格等于T-only冻结资格集合、LOPO移除整个heldpatch；raw结构与stage4 task13封存hash不变。

原dev冻结代码另存`frozen_preconfirmation/qualify_training_dev_v1.py`。Astra在确认前发现未来sphere-invalid会导致矩阵跳行与LOPOmask长度不一致，已改成整系统明确geometry_invalid、无静默删行，并由[test_geometry_path.py](test_geometry_path.py)确定性控制通过。[修订记录](geometry_path_amendment.json)保存前后hash，旧dev无此非法point，未重新求dev数值，旧freeze保留。这是未来失败路径修正，不是用评价结果改拟合。

当前检查点只到开发；已封[训练](real/development_training_freeze.json)、[资格/全部T系统](real/development_qualification_freeze.json)、[评价配置](real_evaluation_freeze.json)、none选择及失败明细。待Astra审阅后再决定固定确认阶段，不追加patch、window、第五方法或108。
