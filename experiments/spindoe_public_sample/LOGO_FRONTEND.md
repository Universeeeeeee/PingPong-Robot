# 无点阵模板的商标纹理前端

已实现：球圆归一化 → 球内暗色印刷纹理角点 → LK前后向跟踪与ID → 可见球面提升 → 相邻时刻有限旋转拟合。

新分支不加载 DotDetector、ref_pattern.csv 或 BayGeoHasher。普通商标没有现成的三维点阵对应，因此不能只替换圆点检测器后继续调用旧的点阵姿态恢复；本分支直接拟合同一轨迹的两次球面位置。旧 SpinDOE 实验保持原样，供对照。

## 本次找到的真实素材

| 来源 | 已获取内容 | 是否可用于测量 |
|---|---|---|
| [Sony官方视频](https://sonyresearch.github.io/ace_public/) | `data/logo/sony_gcs_official.mp4`，93,082,306字节，1920×1080，播放120fps、36.783s，SHA256已校验 | 放大panel是事件可视化，混有后处理旋转轴/白点/文字，镜头还在主动跟踪；未作为干净输入拟合，也未把屏幕转速当真值 |
| [图宾根官方商标研究](https://uni-tuebingen.de/en/132097) | [7格真实商标序列图](https://uni-tuebingen.de/fileadmin/_processed_/0/f/csm_moving_logo_6928586980.webp)，419×59像素；拆成7格做回放 | 可用于定性前端测试；每格实际采集间隔未知，禁止把论文380Hz直接赋给这些展示格 |
| [Pnachtwey 2000fps原帖](https://www.tabletennisdaily.com/forum/topics/how-does-one-post-a-video-on-you-tube.4053/) | 核实作者提供过原始AVI和采集设置 | AVI迁移后404，公开YouTube副本不可用，未获得素材 |

尚未找到可直接下载且同时具有清晰商标、无算法叠加、明确采集时序的完整原始视频。Sony完整文件可作参考演示，不能取代这种测量输入。图宾根回放是由真实静态展示图的7格生成，**不是原始拍摄视频，也没有插值生成中间姿态**。

## 实际测试结果

- 图宾根7格序列：6个相邻窗口中1个通过局部拟合检查，其余因轨迹不足、对应不一致或拟合残差大拒绝。通过窗口估计每格约18.88°，只能作为候选角位移；没有真实转角标签，不能说明准确。没有输出rad/s或rps。
- 合成三维文字球视频：48帧，已知125Hz，带球心移动与尺度变化；47个窗口中45个通过局部检查。真值 `[1.4,-2.2,3.0]` rad/s；中位估计约 `[1.446,-2.121,2.906]` rad/s，每窗口角速度向量误差中位数约0.728rad/s。合成结果不是实拍验证。
- 6项前端测试通过：文字真正映射在三维球面、多轴/正反转、平移和尺度补偿、空白球拒绝、纹理消失后的ID重置、大步不可跟踪的特定例子、非法时间戳。

## 运行

函数接口：`estimate_sequence(frames, timestamps_s, circles)`，输入RGB uint8帧和每帧的球心/半径 `(cx,cy,r)`。可直接连接动捕投影或人工核对后的球圆。

真实序列图回放（不赋真实时钟）：

```sh
.venv/bin/python run_logo_video.py data/logo/tuebingen_strip_replay.mp4 --circles data/logo/tuebingen_circles.json --count 7 --out outputs/logo_real_strip
```

新视频需提供球圆JSON，形如 `{"start_frame":0,"circles":[[cx,cy,r],...]}`；固定球圆可用 `{"circle":[cx,cy,r]}`。只有确认采集间隔后才传 `--capture-fps`。不传时，输出字段为 `rotation_vector_rad_per_video_frame`，绝不读取MP4播放FPS来冒充采集时间。

合成验证回放：

```sh
.venv/bin/python run_logo_video.py data/logo/synthetic_text_ball.mp4 --circles data/logo/synthetic_circles.json --capture-fps 125 --count 48 --out outputs/logo_synthetic
.venv/bin/python -m unittest test_logo_frontend -v
```

输出包括 `results.json`、`tracks_preview.mp4`、`preview.png` 和未叠加预测的 `source_frame.png`。序列图预览每秒播放2格；合成预览每秒20帧，显示速度与测量时钟分开。所有叠加视频都标明输入类别。

## 当前边界

- 这是球内暗色纹理前端，不是OCR或语义商标分割。有反光贴点时仍需商标掩膜/ROI，否则也可能跟踪其它标记。
- 必须有可靠球心/半径；模型采用弱透视，圆参数误差可能被解释为旋转。
- 输出在观测坐标系，相机自身的转动尚未补偿；需要固定观察方向，或另外提供相机姿态来转换到球台坐标。
- 前后向检查及有限旋转拟合只检验局部一致性。`valid`不表示真实角速度已获认证。
- 当前只接受小的**估计**帧间转角，并假设真实运动也足够小。无法排除整圈混叠：每帧恰好转一圈可能被观察为静止，20°门槛不能识别这种情况。
- 商标不可见、轨迹不足或特征分布不足时输出无效，不将无效当无旋。
- 普通帧前端已实现；未经原始G1事件数据验证，不能据此宣称G1上的识别效果。

下一步需要一段未叠加算法图形的真实商标原始序列，连同采集时间、球圆/标定信息，才能做实拍效果和误差评估。现有公开素材不足以作旋转真值系统验收。
