# 2×2局部测量合成gate（2026-11-07）

**clean oracle测量校准通过，tau固定为0.5px。** 合格正控制H/下一1ms短future共74条等权距离记录，最坏P95为0.33544px（reverse patch2、第一段、pol0、edge3、future，5个不同pixel），1.25倍为0.41930px，因此取预定0.5px floor；没有扩大阈值。尚未生成任何新的real模型结果，也未运行CMax或108。

[protocol_synthetic.md](protocol_synthetic.md)明确：每patch三段×两polarity，最多36个局部观测/window；两极性不是独立录制重复。各representation支持相同的M_time/M_spatial四组合；R_all保留事件，R_latest仅本2ms段末按rawpixel取latest，不搜age。M_spatial使用pixel协方差回归，不做混单位3D TLS，允许vn=0；方向不唯一/时间不足保留不可识别。旧timegate只诊断，未用于新模型筛选。

固定stage4六scene，加base的{0,0.5,1}ms正偏移尾迹和fixed-edge同步flicker/同pol时间shuffle两个控制。8scene选格只看前6ms T event数，得到35个合成格，共840条2×2记录（含空/失败），truth只交oracle evaluator。T/H/V先按rawpixel和时间分割，H/短future按(pixel,pol,0.5ms bin)先单位内全event平均，再单位等权，不挑时间对应。

tau只由clean三rot、mixed、puretrans中满足≥12 T unique、H≥4不同pixel或future≥4不同pixel/≥2bin的单edge oracle构成；puretrans当前各2ms段没有足够12个T unique，作为支撑不足控制保留，不虚构它的有效校准分数。因此0.5px只覆盖这些实际可识别的clean局部观测域，不能外推所有慢平移/稀疏格。

clean三旋转/mixed的2×2模型，在≥12 T unique且可评价的短future上，P95中位约0.10–0.13px，R_all和R_latest在干净单crossing下基本一致。尾迹不进tau：oracle短futureP95最高约3.95px；R_all/M_time及M_spatial的可评future P95中位约3.31/1.45px，R_latest约4.66/1.77px，均不能因相对改善称通过0.5px线。它是多层/延迟观测失配，不通过抬tau掩盖。固定flicker同pol短future仅一个半msbin，不满足两bin，保留不可评；它可支持静态constraint，不能当非零motion成功。零scene没有event。

[model_verification.json](model_verification.json)七组正例/单位/退化/同支持/等权测试通过，包括真实vn=0固定line可以被M_spatial辨认、M_time逆时间梯度则不可辨；[verification_synthetic.json](verification_synthetic.json)独立复算全部支持和计分、oracle隐式边缘有限差分、冻结来源hash均通过。oracle normal最大差1.03e−10，speed最大差1.31e−6px/s；单位P95差0。原始结构检查点和stage4 task13封存未改。

目前模型可用只表示代数可辨，不等于局部非零motion资格。后续training资格只能用T支持/时间跨度/32个pixelcluster bootstrap等训练量；H/shortfuture/19 permutation只评价运动证据，**不能筛掉local约束后再拟合共同omega**。stationary与nonzero分开，预测位移低于tau不可辨。完整静态/permutation/global/leavepatchout规则待Astra真实协议冻结，未提前应用或用real选参数。

交付：[oracle校准](oracle_calibration.json)、[全部oracle记录](oracle_observations.json)、[全部2×2模型](model_observations.json)、[输入/格选择](input_manifest.json)、两个小型控制`inputs/`及独立truth。H支撑的最终“不同pixel”口径按Astra在real之前澄清，从已保存单位键复核，tau保持0.5；原校准摘要另存prereview版本。到此交合成检查点，等待真实阶段gate。
