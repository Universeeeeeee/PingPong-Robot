# B2：不规则商标候选与少量圆样辅助观测

旧9个人工种子仍不能可靠指认商标；将circle失败或薄弧叫logo会误判。全片固定survey后，在11.75s发现具有连续非单薄弧细节的不规则商标候选L。新片段为post-survey探索性增补，符合用户形态提示，尚未排除圆点斜视、重叠和ROI截断，不能称物理商标身份已确认。

L在11.750–11.800s的5ms共10/10块、2ms共25/25块通过B0原匹配gate（均含1人工seed），是局部形态/关联延续证据，不是商标识别准确率。然而减一次初始偏置后L的2/5ms相对位移差仍最大9.520px，N1最大7.977px，尚未达到稳定几何位置观测。辅助点仍独立跟踪，既不改变也不回填L；联合门槛只改变整组观测覆盖，没有证明中心、姿态或轴准确性提高。

## 人工提示、原始证据与选择过程

用户提示“不规则外观=商标、圆形=贴点”用于人工候选分层；并非以任意形态阈值训练分类器。旧主样本固定不替换。先审[参考全ROI](reference_context.png)和card_*局部100ms演化卡：上排5ms帧0/5/10/15/19，下排最近2ms原始窗，mid差±0.5ms。失效时明确固定seed审阅，绝不补成同一特征轨迹。

父代理按全片0.25+0.5k的35个固定5ms窗[定时survey](whole_recording_survey_1.png)（另whole_recording_survey_2.png）检索，随后审[三个连续窗细节](irregular_survey_details.png)，再提出11.70–11.80s探索分支。新种子L(582,212)、较远A1(533,192)、A2(523,238)先于新跟踪冻结；之后父代理仅看原始seed图、尚未见跟踪结果时提出near N1(588,190)、N2(576,237)，单独追加near_protocol.json，原三条输出逐字保持。新数据选择不是独立测试，不能据此估计泛化率。

先看[完整5ms原始序列](exploratory_raw_5ms.png)、[完整2ms对照](exploratory_raw_2ms.png)、[跟踪前原始候选卡](exploratory_seed_review.png)，再看[较远对照跟踪](exploratory_tracking_review.png)与[near跟踪](near_tracking_review.png)。L框有多方向细节，并非仅单一细弧，但上方圆弧/右侧相邻结构仍可能混入。

## 候选证据卡：人工推断与几何gate分开

旧样本下表为100ms的20个5ms块；新分支为向前50ms的10块。circle/ellipse数量仅是原B0局部模型可拟合gate，不能验证语义；连通分量是事件碎片，不是Marker数；边框占比描述27px截断风险，不是物理遮挡真值。

|片段/F|人工候选层级|有效块/总块|circle/ellipse gate有效|人工依据/替代解释|
|---|---|---:|---:|---|
|early_turn/F1|round_candidate|20/20|3/6|双极性圆样边缘可见；27px框截断/内部填充使圆模型支持随时间变化|
|early_turn/F2|ambiguous|20/20|0/1|混合多弧、初始主要位于框外；不可由不圆或拟合失败推商标|
|early_turn/F3|round_partial_candidate|13/20|5/8|起始圆/椭圆边缘，随后变窄和离开ROI；丢失后不作同一结构解释|
|middle/F1|round_candidate|20/20|11/16|原始圆样边缘连续，部分窗口圆/椭圆模型支持|
|middle/F2|round_candidate|20/20|3/3|圆样双极性结构，锚点框常截去一侧；模型少不代表logo|
|middle/F3|ambiguous|19/20|1/5|薄弧/稀疏极性结构，不能排除斜视圆点或单一边缘|
|late/F1|round_candidate|20/20|18/18|分散的双极性圆/椭圆边缘，多数局部模型gate有效|
|late/F2|round_candidate|20/20|18/19|分散圆/椭圆边缘，多数模型gate有效|
|late/F3|round_partial_candidate|20/20|16/19|近下方截断的圆/窄弧候选，边界/遮挡风险保留|
|post_survey_11_75/L|irregular_logo_candidate_UNCONFIRMED|10/10|0/0|全片survey后选取；多个连续原始块有非单薄弧的多方向细节，符合用户形态提示；邻近圆弧/重叠尚未排除|
|post_survey_11_75/A1|round_candidate|8/10|3/7|与L约53px的较远圆弧对照，不称附近|
|post_survey_11_75/A2|round_partial_candidate|6/10|0/2|与L约64px的较远窄弧/下方截断对照|
|post_survey_11_75/N1|round_candidate|10/10|0/0|与L约23px的附近圆/椭圆候选，27px支持可混入L|
|post_survey_11_75/N2|round_partial_candidate|10/10|0/0|与L约26px的附近窄弧候选，可能共享局部事件|

所有原始形态指标及失败保留在morphology.csv；候选层级与人工依据在candidate_evidence.csv。早F2混合/截断、middle F3薄弧均保持ambiguous；middle F2圆样边缘常被27px框截掉一侧，说明模型失败不能推出logo。新L的圆/椭圆gate均失败也不是判定商标的依据。

## 新分支固定门槛前向50ms结果

仍为27×27双极性不同像素occupancy、sigma0.8、±12px整数平移、≥12像素、cosine≥0.35、次峰（排除3px邻域）分差≥0.025、FB≤2px、不触搜索边。参考块5ms为[11.750,11.755)，2ms为[11.750,11.752)，不是同一瞬间；一次拒绝后不接回。没有固定贴点布局、共同平移或平面假设。

|特征|T ms|有效/总块|circle/ellipse gate有效|首次拒绝|
|---|---:|---:|---:|---|
|L|5|10/10|0/0|none|
|L|2|25/25|0/0|none|
|A1|5|8/10|3/7|rejected_ambiguous_peak+forward_backward|
|A1|2|12/25|4/11|rejected_ambiguous_peak+forward_backward|
|A2|5|6/10|0/2|rejected_forward_backward|
|A2|2|4/25|1/3|rejected_ambiguous_peak|
|N1|5|10/10|0/0|none|
|N1|2|25/25|5/14|none|
|N2|5|10/10|0/0|none|
|N2|2|7/25|1/2|rejected_ambiguous_peak|

L和N1的2ms位置只在相邻有效块间线性插值到5ms mid；相对位移差仅减最接近名义参考时刻的一个固定初始偏置，绝不逐帧重对齐、不跨失效插值。下面为中位/P95/最大欧氏差，单位px；初始零相对差由定义产生，不是额外真值。

|特征|配对数|原始锚点差 中位/P95/最大|减单个初始偏置后的位移差 中位/P95/最大|
|---|---:|---|---|
|L|10|6.951/9.975/10.177|6.375/9.312/9.520|
|N1|10|5.770/7.462/7.603|5.817/7.809/7.977|

整数位移量化、模板更新漂移、窗口内容和ROI截断均可贡献差异，不能解释为商标/贴点中心误差；完整位置/初始参考时间见cross_duration_positions.csv，变化图见[cross_duration_differences.png](cross_duration_differences.png)。

## 公平比较：单候选与联合可用性

每组使用完全相同主体轨迹，不平均不同点的2D位置/速度，不要求投影距离恒定，不用未验证辅助补主体。position联合gate要求指定所有轨迹当帧可用；aux-model联合gate另外要求每个辅助至少一种circle/ellipse gate有效，主体L不强套圆模型。单主体的辅助模型列为NA，不代表L几何中心已有效。

|分支|T|组|主体可用|联合position可用/总块|联合aux-model可用|主体有效时因辅助失去的覆盖块|
|---|---:|---|---:|---:|---:|---:|
|new_distant_control|5|L|10|10/10|NA|0|
|new_distant_control|5|L+A1|10|8/10|7|2|
|new_distant_control|5|L+A1+A2|10|6/10|2|4|
|new_distant_control|2|L|25|25/25|NA|0|
|new_distant_control|2|L+A1|25|12/25|11|13|
|new_distant_control|2|L+A1+A2|25|4/25|3|21|
|new_NEAR|5|L|10|10/10|NA|0|
|new_NEAR|5|L+N1|10|10/10|0|0|
|new_NEAR|5|L+N1+N2|10|10/10|0|0|
|new_NEAR|2|L|25|25/25|NA|0|
|new_NEAR|2|L+N1|25|25/25|14|0|
|new_NEAR|2|L+N1+N2|25|7/25|0|18|

旧主样本条件组压缩如下（5ms向前50ms；2ms逐块结果留CSV），主体均为ambiguous假定不规则，不能称商标真值。

|旧条件组|单主体position|+1辅助position / aux-model|+2辅助position / aux-model|
|---|---:|---:|---:|
|old_ambiguous_early|10/10|10 / 6|3 / 0|
|old_ambiguous_middle|9/10|9 / 7|9 / 2|

旧early/middle两组主体仍是ambiguous假定不规则候选，不能称商标对照真值。新near组更贴近“商标附近点”的构想；A1/A2距L约53/64px，仅为较远分散round对照，N1/N2约23/26px。辅助轨迹提供其自身事件/形态延续和局部模型质量证据，不自动降低L匹配歧义。联合gate是逐帧质量标记，可因模型gate波动重新有效；匹配真正丢失的ID始终不复活。

5ms near N1/N2虽均10/10块可跟踪，但其局部圆/椭圆模型均0/10有效；同时有两条辅助轨迹不等于增加两个可信几何中心。patch模板更新可能沿边缘游移，尤其N2锚点不等于可见圆/弧的中心。2ms辅助模型有部分支持，但N2较早因歧义丢失，联合覆盖相应下降。

## near支持重叠：不能称独立几何约束

|辅助|T|同时有效块|最小锚点距离px|最大27px patch IoU|最多共享原始事件/不同像素|最多共享L事件比例|
|---|---:|---:|---:|---:|---:|---:|
|N1|5|10|22.804|0.078|46/41|17.6%|
|N1|2|25|20.248|0.123|31/31|36.9%|
|N2|5|10|24.331|0.050|17/16|6.5%|
|N2|2|7|23.195|0.070|10/10|10.8%|
|A1|5|8|49.980|0.000|0/0|0.0%|
|A1|2|12|48.270|0.000|0/0|0.0%|
|A2|5|6|62.177|0.000|0/0|0.0%|
|A2|2|4|63.071|0.000|0/0|0.0%|

共享支持直接由同一原始窗的事件索引/像素交集计数，说明邻近patch不能直接当独立轴约束。IoU为零也不证明不同物理身份或独立统计噪声。本轮没有由多轨迹数量推导三维约束或识别准确性。

## 仍不能确定与下一步

用户形态提示在新L原始事件上有可审阅落实：稳定不规则候选可按原门槛短时跟踪；旧9样本不能可靠确认商标这一负结果保留。尚未确认L确为商标而非混合/投影结构，也未确认patch锚点是同一物理球面点或几何中心。

目前不能宣称辅助改善主体匹配：方法仍是独立跟踪，联合筛选会损失覆盖，附近patch共享事件且5ms局部圆中心支持不足。辅助轨迹没有用于修正L，因此本轮也不能判断经验证观测模型下的辅助融合能否提升精度。需要人工跨帧确认商标/圆点边界与重叠、真实遮挡/出入视野及不同布局记录，才可评价身份与失效检测；下一步先分离整数位移量化、递推模板累积漂移与ROI混入的贡献：固定参考外观模板对递推模板、亚像素对整数位移、不同ROI支持各作独立对照，不同时优化。再检查L内部可重复角点/笔画交汇及辅助Marker独立边界，验证何者对应同一物理位置；不急继续加点，不能将整个patch锚点直接提升到球面或送PnP。需要经验证的观测模型才能利用辅助构造刚体或姿态约束。PnP、三维旋转与轴准确性均未启动。

## 复现与核验

旧630行、90个审阅panel、210个绝对时间窗保留；363个有效形态行由原始xy集合和独立8邻域flood-fill/边框坐标核对。新175行全保留；适用于2/5ms长度的前向循环实现与原B0的315行完全相同；15个新接受/拒绝匹配25×25分数图与次峰直接内积复算通过。L/A1/A2原三条结果在near追加后hash未变，全部10路径丢失后未重连，联合条件不改主体。语义身份没有被这些实现检查验证。

```sh
g1_spin/.venv/bin/python g1_spin/experiments/logo_local_20261102/review.py
g1_spin/.venv/bin/python g1_spin/experiments/logo_local_20261102/validate_review.py
g1_spin/.venv/bin/python g1_spin/experiments/logo_local_20261102/exploratory.py
g1_spin/.venv/bin/python g1_spin/experiments/logo_local_20261102/exploratory.py --track
g1_spin/.venv/bin/python g1_spin/experiments/logo_local_20261102/nearby.py
g1_spin/.venv/bin/python g1_spin/experiments/logo_local_20261102/compare.py
```

父代理survey复现脚本whole_recording_survey.py与survey_followup.py另存本目录。协议为protocol.md、exploratory_protocol.json、near_protocol.json；原始来源part hash见review_metadata，脚本/配置hash见metadata与validation及各分支metadata。未改原始数据/B0/B1/旋转，未训练、联网、装依赖或操作硬件。
