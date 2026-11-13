# 局部运动测量合成检查点（2026-11-07）

最终计数口径（Astra在新real之前确认）：每patch三训练段×两polarity分别处理，每窗最多36个局部时空观测，七窗最多252个位置/时间/极性观测×四个2×2组合；空与失败保留，双极性不当独立录制重复。下文“每格至多三个训练段”并非只保留一种极性。

仅合成metric/oracle校准及2×2实现核验；禁止新real结果、CMax、108回归或自动改门限。raw结构检查点保持封存。Astra high-level方案、Sol执行；实录完整协议及permutation/静态决策终稿尚待下发，不将此中间实现冒称最终实录资格。

## 固定表示和模型

每个12×12 raw sensor格至多三个不重叠训练段[0,2)、[2,4)、[4,6)ms。先按rawpixel分T/H：`(x+2*y)%5==0`为H，其余T；V为[6,10)ms所有pixel。时间与pixel分割先于表示，T/H/V不混置乱。R_all保留该段同polarity原events；R_latest在该2ms段末因果地按rawpixel保留最近T event（并列按最大物理ID）。同一representation的M_time和M_spatial使用完全相同支持，不改邻域、按残差删除event或补失败格。

M_time最小化时间残差t=g·(x−x0)+c；转n=g/||g||、vn=1/||g||，normal与speed方向成对处理。M_spatial以Cxx−CxtCtx/Ctt的最小特征向量作为单位n，vn=nᵀCxt/Ctt，直接最小化pixel法向残差；不做任意px/ms混单位3D TLS。若Ctt不足、空间支持不足或最小方向不唯一，记录不可识别。此合成实现检查使用Ctt>1e−12s²、eigengap>1e−8×max(1,trace(Cxx))的数值可辨阈值，仅防代数退化，不作real拟合门限优化。符号以n首个非零分量为正，vn/intercept同步翻转。

每obs至少12个T unique rawpixel才能进入后续局部模型资格；H≥4个不同pixel，future≥4个不同pixel且≥2个0.5ms bin。计分仍按pixel/pol/bin等权，而非用同pixel多个bin替代空间支撑。旧time RMS/std门限仅报告，不作为新空间模型接受门限。pixelcluster bootstrap32、normal散布、速度区间/偏差和19固定timepermutation将在决策终稿前只作固定实现准备，不得据实录调参。

## 固定合成输入与观测控制

沿用stage4 base/reverse/second/translation/mixed/zero完整整数event；不重新生成或修改其输入。额外只有：

1. base同pixel尾迹控制：保留原event，并追加+0.5ms、+1ms、仍在10ms窗内的时间偏移event，pixel/polarity保持。truth只在evaluator记录base来源/root及偏移，runner不能用偏移补偿。
2. fixed-edge synchronous flicker/time-shuffle控制：一条固定竖直edge周围整数pixel在固定0.5ms周期上同步响应，随后按固定seed7将时间在相同polarity中置乱；不引入真实平移。该控制用于检查静态constraint和伪运动，不能算非零运动正例。

clean三rot、translation、mixed为具有明确局部几何motion的正校准来源。tail和static flicker/shuffle作为观测失配/负控制单列，不通过抬高tau掩盖时间尾迹或同步触发。zero无event保留不可测，不算测准zero。

每clean/尾迹合成scene按前6ms T event数，在12×12格中固定count降序/坐标排序，≥20events、格中心距≥16px、最多6格，之后不按residual或truth挑格。每格完整event保留物理ID。oracle只有evaluator可按材料edge标签分组：在obs段中心，使用生成器材料normal、已知omega和投影球面点推导真实瞬时pixel n/vn，只从该edge T events估intercept，评价全部同edge H/短future。多个edge混杂保留各edge结果，不让runner使用truth分组。

本合成oracle使用真值分组只校准距离测量分辨率，不作为实录patch合格前提或runner选择。runner则对每个patch/polarity全部T拟合，不读edgeID/omega。有限曲率/变速导致局部线性近似误差如实留在oracle中。

## H与短future度量

空间法向距离为`abs(n·(xy−x0)−vn*(t−t0)−intercept)`，不按预测残差选择event。H取当前2ms训练段的H rawpixel。短future只取下一1ms：对[0,2)段评价[2,3) H像素；[2,4)段评价[4,5) H像素；[4,6)段评价[6,7)全部V。每个(pixel,polarity,固定0.5ms bin)先平均该单位全部event距离，再单位等权统计/P95/coverage；不挑最近时刻或对应。未来不足4单位/2bins明确不可评，不删掉坏单位。

在支持条件满足的clean正控制oracle H和短future距离中，tau=max(0.5px,1.25×最大P95)。若tau≥3px，说明此固定尺度的测量校准失败，报告具体scene/patch/edge/segment和分辨率，不增阈值、不扩大格或换控制。tail/staticoracle与events估模型均报告，但不混入clean正控制tau。校准、实现测试和manifest完成后交Astra gate，再等最终实录协议。所有原始来源和旧结构诊断hash保持。
