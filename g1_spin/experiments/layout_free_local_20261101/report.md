# B1：布局无模板局部关联压力测试

结论：B0本来就不需要三维布局模板，但仍依赖局部外观模板。原三段各100ms数据中，本轮只评估参考frame10向未来frame19的50ms。三例目标patch事件移除均在首个干预块拒绝且不复活；+8px同帧复制混合未触发预声明复制锚点吸引，但造成最大1.414px相对基线偏离，其中一例在副本撤去后歧义拒绝。这只是固定输入压力证据，尚不能称可靠拒绝一般干扰、商标识别或跨真实布局泛化。

没有人工商标/Logo标签，没有真实换布局记录。所有F仍是局部特征名；基线锚点不是物理中心或身份真值。复制相同事件patch也不等同真实新增圆贴点，不能据此否定不规则商标可能提供的额外区分力。

## 固定协议与评估范围

预声明见[protocol.md](protocol.md)。保留B0 9种子、5ms、27px双极性binary occupancy外观patch、sigma0.8、±12px整数搜索、12不同像素、score0.35、次峰差0.025、FB2px与不触搜索边的全部门槛。仅前向时间更新外观模板；一次拒绝即停，原信号恢复也不在新位置复用旧F。未增加固定参考模板或亚像素改进。

9条基线90行的状态、锚点、最佳/次佳分数、FB等匹配字段逐位完全复现B0的5ms前向结果；每F单独重复调用亦相同。它们原本就独立，删除其他seed条目并没有删除邻居事件，也不是换布局实验。前半段仅作原始背景，不拿反向传播补可见性。

人工参考块为early [1.750,1.755)、middle [6.350,6.355)、late [11.350,11.355)，各向前到1.800/6.400/11.400，共10块含seed、9次自动匹配。正常匹配只用来源和当前事件；FB也只回查来源块，不读取未来块。

## 输入压力与oracle边界

固定代表为三个F1，frame13–15即名义中点后15–30ms施加压力。移除/复制位置由B0完整基线轨迹构造，属于非因果oracle输入压力，不用于正常算法初始化/预测/控制。所有原始输入只读，修改后的副本归档为*_events.npz，可逐条核查added/removed/clone_source。

- `occlusion`文件名仅为简写：实际移除baseline锚点±13px方框内事件，并非遮挡一个已识别物理Marker。方框外原始事件逐条保留，邻居重叠部分可能随方框一起被移除。
- `clone_plus8`复制同帧该方框事件，x+8px，y/t/p及其他字段不改；原事件全保留。两个patch可混合，并非两个独立Marker投影。重复时刻/同像素按B0 binary occupancy归一化，不用raw新增计数冒充有效干扰。
- `outside_context_delete`是有限实现隔离检查：按baseline每帧当前/前一锚点±30px并集保护，删除外部事件。局部外观模板之外，±12px搜索与sigma0.8/truncate4的3px核使源/目标实际上下文至少±28px；30px保守保护并覆盖两相邻步。匹配图像边界x450..670/y120..280不变。此oracle检查不是布局泛化成绩。

## 干预结果：所有尝试与停止保留

|片段|输入|有效/10（含seed）|首拒绝帧/原因|接受锚点相对基线最大偏离px|距复制锚点最小值px|停止后帧数|
|---|---|---:|---|---:|---:|---:|
|early_turn|occlusion|3|f13 rejected_forward_backward|0.000|NA|6|
|early_turn|clone_plus8|10|none|1.414|7.000|0|
|early_turn|outside_context_delete|10|none|0.000|NA|0|
|middle|occlusion|3|f13 rejected_insufficient_pixels+score+ambiguous_peak+forward_backward|0.000|NA|6|
|middle|clone_plus8|6|f16 rejected_ambiguous_peak|1.414|7.071|3|
|middle|outside_context_delete|10|none|0.000|NA|0|
|late|occlusion|3|f13 rejected_insufficient_pixels|0.000|NA|6|
|late|clone_plus8|10|none|1.414|8.000|0|
|late|outside_context_delete|10|none|0.000|NA|0|

三个patch移除例均在frame13首个受干预块拒绝；frame16原始事件恢复后仍不自动复活。这支持当前三例的拒绝/不接回行为，不能将3/3当部署成功率，亦不能确认单个物理Marker被遮挡。图中邻近结构仍存在，neighbor_patch_counts.csv记录邻居baseline patch计数：并非把整球或整个搜索区抹空。

复制混合时三个F1在frame13–15仍接受，最大偏离1.414px；早/晚偏离在复制撤去后持续，中段frame16以歧义拒绝。逐帧模板更新可能继承混合造成的偏置，但没有独立身份标签，不能据此判真实identity switch或绝对定位误差。所有接受位置均未满足“距复制锚点≤3px、距基线≥5px且更近至少3px”的预声明copy_anchor_attraction条件；未发生该几何指标不证明一般身份安全。

以下为首受干预块frame13及副本撤去frame16的实际匹配证据；全部10块、拒绝提议/次峰坐标保存在observations.csv。

|片段|输入|帧|状态|提议x,y|次佳x,y|score/margin|FB px|
|---|---|---:|---|---|---|---|---:|
|early_turn|occlusion|13|rejected_forward_backward|555,207|555,211|0.426/0.0332|2.236|
|early_turn|clone_plus8|13|accepted|553,208|559,208|0.819/0.0918|1.000|
|early_turn|clone_plus8|16|accepted|541,209|540,213|0.802/0.1116|0.000|
|middle|occlusion|13|rejected_insufficient_pixels+score+ambiguous_peak+forward_backward|546,228|543,219|0.320/0.0020|2.236|
|middle|clone_plus8|13|accepted|534,224|542,225|0.700/0.0685|0.000|
|middle|clone_plus8|16|rejected_ambiguous_peak|529,226|522,227|0.595/0.0225|1.414|
|late|occlusion|13|rejected_insufficient_pixels|533,225|534,229|0.398/0.0821|1.000|
|late|clone_plus8|13|accepted|529,221|537,221|0.678/0.0790|0.000|
|late|clone_plus8|16|accepted|525,220|525,224|0.637/0.0914|0.000|

## 原始新增事件与归一化后实际变化

表中每例为3个活跃干预块的总raw添加/删除事件，以及每块改变的极性×像素occupancy数范围；ON和OFF在同像素分别计通道。这是binary输入实际改变，不是Marker数量。

|片段|输入|raw新增/删除总计|每块改变occupancy范围|
|---|---|---:|---:|
|early_turn|occlusion|0/791|225–264|
|early_turn|clone_plus8|791/0|201–227|
|middle|occlusion|0/418|111–151|
|middle|clone_plus8|418/0|110–147|
|late|occlusion|0/323|96–116|
|late|clone_plus8|323/0|94–114|

`outside_context_delete`三个F1全部输出与基线相同；源/目标正反向所访问的108个blur数组切片逐位相同。其证明这些特定输入路径对实际上下文外变化不敏感，是预期的局部实现隔离属性。不能推导邻近结构（能进入搜索/blur上下文）任意变化也不影响关联，更不能证明真实换布局泛化。

## 阶段出口与所缺数据

当前可以保留：无需三维布局模板的人工种子局部关联；三例强patch证据移除能及时拒绝并维持丢失；上下文外更改的有限隔离。不能保留的宣称：已检测商标、已识别物理Marker、真实换布局鲁棒、一般可靠拒绝、三维轴准确。局部外观模板仍存在，逐帧更新产生偏离的机制尚未排除。

下一步需要人工确认商标/特征身份的连续原始记录、同球或不同球真实改变邻近贴点布局的重复记录，以及短时真实遮挡/出视野/邻居交叉并有人工身份核对的片段。只有这些数据才能区分不规则商标的区分力、真实布局变化与合成叠加。本轮不新增算法，固定参考模板对逐帧更新可作为独立下一轮消融；先保留本轮全部失败与偏离作为对照。

## 复现与独立检查

基线90行、干预90行（3片段×3输入×10帧）完整唯一；90个归档事件payload直接用原SDK绝对时间mask和独立方框定义核对，clone的t/p及其他字段不变，原事件全保留；occupancy差用独立(x,y,polarity)集合核对。12个关键匹配的25×25位移分数图和次峰以直接内积独立复算，均一致。108个源/目标blur访问切片逐位相同；全部9条压力路径核对拒绝后不重连。旧B0目录文件hash未变，原始输入manifest/part、协议、脚本、副本hash在metadata.json；验证见validation.json。

```sh
g1_spin/.venv/bin/python g1_spin/experiments/layout_free_local_20261101/run.py
g1_spin/.venv/bin/python g1_spin/experiments/layout_free_local_20261101/summarize.py
```

图：baseline_review.png原始参考；occlusion_review.png目标patch移除；clone_plus8_review.png合成复制混合。灰+为oracle基线锚点，洋红×为接受锚点，红×为拒绝提议；失去后不补点。表格包含baseline、observations、interventions、summary、context_checks、neighbor_patch_counts、independent_feature_checks。仅写此新目录，无PnP/旋转/训练/联网/硬件，不安装依赖。
