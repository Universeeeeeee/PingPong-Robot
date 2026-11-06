# 阶段5最终报告（2026-11-07）：时间面结构冲突已定位，局部运动仍未通过确认

**本轮固定原始诊断、合成测量校准、2×2开发和三窗确认均已完成，selection保持none。** 旧F0的大部分失败可由重复pixel/polarity的多值时间结构直接解释；pixel法向拟合改善了部分训练方向稳定性，但没有局部运动同时通过留出距离、静态对照、置乱与可辨分辨率。因此本轮未获得可验证实录角速度，按停止条件没有CMax或108扩展。

本报告使用**修正完整raw输入后的确认结果**。Astra终审发现首版确认从旧common_ids填充12px格，边界格被截断，与开发的完整官方解码格定义不一致；首版确认、旧报告及其hash均保留`local_motion/prereview_clipped/`和`report_prereview_clipped_confirmation.md`，不作为最终确认结论。保持原格、几何、tau、bootstrap/permutation种子和规则，独立重建完整raw物理ID后仅重跑run3/run5；run4全部格输入及144条训练/144条评价精确一致复用。

目录保留`drill_phase5_20261106`启动标识，实际完成日期2026-11-07。方案与各检查点由GPT-6 Astra xhigh审阅，GPT-6.1 Sol high执行。旧源码、RAW和阶段1–4证据保持，不根据确认结果改模型或门限；完整RAW不上传。

## 1. 原始局部结构：先解释旧门限为什么不可能通过

按前6ms事件数选固定12×12px格，不看残差或omega：run1 1.495s、run2 2.095s、run5 4.095s分别6/6/3格，静止1.525s为0格。15格完整原始10ms event的极性分离x/y/t、半ms切片、pixel计数/raster、first/median/latest时间面与同pixel间隔均保留，见[原始结构报告](report_raw_structure.md)。后4ms只可视观察，没有用于训练/选格。

旧gate RMS≤max(20us,0.15×std_t)，非floor分支相当于R²≥0.9775。把同pixel同pol的时间variance分成within/between，任何单值t(x,y,pol)的最小RMS也至少是sqrt(within)：旧三个F0开发窗≥10点邻域中，256/269、173/186、71/95的这个下界已超过gate，中位下界/gate为3.9114/3.3404/3.5290。执行者独立重建与Astra低层审阅吻合，说明这些**全部重复event**不可能由更充分的单值时间面优化挽救。

这只解释F0，不解释去重复F1仍失败。同pixel重复主要间隔0.74–3.22ms，所选格≤50us间隔比例0；不能简单归为密集短burst。多个时间层、异极性空间分离与纹理/曲率/变速等解释仍有歧义。球内外事件率无统一强同步，但球外含杆/手/背景，不能据此排除照明，也不能宣布硬件无flow。

## 2. 先用oracle校准测量，再固定2×2

沿用stage4 clean三旋转、纯平移、mixed、zero，再固定尾迹和flicker/time-shuffle两个观测控制。clean可辨单edge oracle的H/下一1ms future共74条等权记录，最大P95=0.33544px；预定tau=max(0.5,1.25×最大P95)因此为**0.5px**。该纯平移控制在本次选格和2ms支持规则下没有任何合格oracle记录进入τ校准；74条实际仅来自base/reverse/second/mixed；尾迹、多层和静态控制不用于抬tau。具体范围和独立有限差分核验见[合成gate](local_motion/report_synthetic.md)。该校准不保证所有实录薄层都适用，更不是物理精度真值。

四组合为R_all/R_latest×M_time/M_spatial：同representation支持，M_time为时间回归；M_spatial用Cxx−CxtCtx/Ctt最小特征向量求单位normal、vn=nᵀCxt/Ctt，最小化pixel法向残差，不作混px/ms的3D TLS。R_latest仅本2ms段末取rawpixel latest，不搜age；vn=0允许stationary约束，近零不除vn。旧timegate只诊断，不再机械放宽。

每patch三段×两polarity，最多36局部obs/window；训练前先分T/H/V，trainer只打开独立T文件。至少12 T unique、≥24/32有效pixelcluster bootstrap、normal角P95≤15°等T规则决定训练资格。candidate_nonzero还需|vn|>max(50px/s,2×bootstrap95偏差)；static-compatible单列，S系统保留实际vn，不按标签赋0。

H/future所有(pixel,polarity,halfmsbin)等权，H≥4不同pixel、下一1ms future≥4且2bin；moving两边P95≤0.5px，并相对同T最佳static line的mean距离改善≥20%，真实min gain须严格优于至少15个可评的19次各split内部时间置乱。位移低于tau标不可辨。H/future只评价，**不筛约束再拟共同omega**。全部T合格系统先hash保存，再读取评价数据。

## 3. 开发：训练方向改善，但验证证据没有跟上

开发3窗+静止控制共15格、360条四组合。run1 Rall/Rlatest的M_spatial分别有14/13个T合格obs、5个patch；run2只有8/6个，平移4/3个、2patch。M_time更少。所有组合的留出支持非零为0、stationary为0，[开发报告](local_motion/report_development.md)和[分层失败](local_motion/real/development_failure_attribution.json)同时给边际通过及串行首败数量，避免用一句allzero掩盖原因。

run1空间模型H P95中位2.16/2.62px，短future3.59/3.70px，仍明显高于0.5px；run2约2.88/2.84和3.37/3.29px。部分方向训练稳定或相对static有gain，但不是完整运动证据。两个run1空间组合虽然满足T-only全局数量/秩要求，全部约束LS残差约498/507px/s，仍只是诊断。

冻结LOPO同heldpatch源模板和同V目标比较omega/2Dtranslation/zero，没有按预测选目标。各解释距离优劣混合，且该旧末1ms模板/3px双向Chamfer未在阶段5重新校准，仅作诊断，不能自动视作global pass或运动证伪。开发组合资格未达，none选择在打开确认之前已封存。

## 4. 固定三窗确认：未出现跨patch支持，不进入CMax

确认run3 1.495s、run4反转2.495s、run5平移5.095s，按相同计数规则选6/6/4格，384条四组合；失败不补格，双极性不作独立录制重复。先T模型/bootstrap/资格/系统hash，后H/future/19置乱。

|确认窗|Rall/time T合格|Rall/spatial T合格|Rlatest/time T合格|Rlatest/spatial T合格|最终非零 / stationary|
|---|---:|---:|---:|---:|---:|
|run3 1.495s|5（3patch）|13（5patch）|5（3patch）|10（5patch）|0 / 0|
|run4 2.495s反转|3（3patch）|8（5patch）|3（2patch）|6（5patch）|0 / 0|
|run5 5.095s平移|1（1patch）|7（3patch）|1（1patch）|6（2patch）|0 / 0|

确认run3的Rall/spatial是唯一T-only全局充分系统，但没有≥3patch的H/future支持nonzero；其余系统不足数量/支撑。修正后run3空间模型H P95中位2.87/2.81px、future2.69/3.12px；反转约2.63/3.15和3.25/3.64px。平移支撑不足/速度不可辨，不能称负例成功。少数观测单项gain或置乱通过，都未同时满足T分类和距离线。

[confirmation_failure_attribution.json](local_motion/real/confirmation_failure_attribution.json)按训练分类、H/future支撑、P95、gain、permutation与位移分别给边际和首败；[confirmation_evaluation.json](local_motion/real/confirmation_evaluation.json)保留全部记录。由于必要条件失败，按预先停止规则**不运行CMax、不做108、不换组合、不增加patch或window**。最终选择保持开发none，并未把确认数据用于重新选择。

## 5. 核验、修订与解释边界

[完整raw输入纠正](local_motion/confirmation_input_correction.json)独立重现run3 patch2/3/4事件187→249、117→259、106→233，run5 patch0/1/3为94→100、107→117、51→79；新增ID及各T/H/V计数逐一保存，其余格不变。旧确认归档，run4输入/结果精确复用；修正不改变开发选择。

[确认增量验证](local_motion/real/verification_confirmation.json)已通过：3窗、384确认+360开发=744个真实2×2模型记录；同representation两fit同支持，T/H/V物理ID互斥；32bootstrap和19置乱计数完整，752次确认等权mean/P95独立复算一致；global/LOPO约束均来自已封T资格，none及旧开发输出不变。此前[原始](verification.json)、[合成](local_motion/verification_synthetic.json)、[开发](local_motion/real/verification_development.json)核验均有独立记录。

确认前已审的sphere-invalid未来分支改为整系统显式geometry_invalid，旧代码/freeze保留，dev无非法point故未重算。确认首次评价在旧hash guard处停止，尚未读H/V，随后只改用确认专属post-T evaluation freeze；阈值/模型/seed和T数值均未改，见[hash路径修订](local_motion/confirmation_evaluator_hash_amendment.json)。这是可追溯配置边界修正，没有用结果调参。

本轮支持：重复event与严格单值时间面门限存在结构冲突；pixel法向模型改善部分T方向稳定性。当前不支持：把这种训练改善称为可验证局部运动、真实旋转恢复或Sony官方复现成功。多边缘/多层、较长多阈值响应、曲率/变速和采样/噪声的具体贡献仍不可唯一识别，真实绝对轴/RPM、zero置信界和在线性能仍无独立真值。

全部既定原始观察、合成校准、开发及三窗确认至此结束。最终源码/完整失败统计/小型裁剪和hash清单保留用于PR4紧凑更新；完整RAW和超大文件不上传。不自动展开下一阶段。
