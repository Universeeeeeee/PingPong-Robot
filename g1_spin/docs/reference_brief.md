# DVSLume-G1 固定球旋转检测：阶段一完整参考资料

我重新核对了 Sony 的论文、官方 SDK 文档和相关开源项目。阶段一的核心资料可以分为四类：算法论文、官方相机 SDK、可复用代码、数学优化与实验验证资料。

这次有一个重要发现：Sony 官方其实公开过一份包含 `ContrastMaximization`、球面事件映射和旋转补偿的 Python 伪代码。 它不是完整可运行的 s-CMax 项目，但对 Codex 复现很有价值，应加入之前的任务书。

## 一、核心算法论文（优先阅读）

1\. [Event-based Gaze Control System for Accurate Real-time Spin Estimation in Professional Ball Games](https://arxiv.org/html/2606.26780v2?utm_source=chatgpt.com)

核心主线

Hu et al. · ECCV 2026 · Sony AI

阶段一最重要的论文。直接对应离线光流初始化、球面 s-CMax、优化策略和固定球电机验证。

重点章节： §3.3.1、§3.3.2、§4.2。

2\. [Table Tennis Ball Spin Estimation with an Event Camera](https://openaccess.thecvf.com/content/CVPR2024W/CVsports/html/Gossard_Table_Tennis_Ball_Spin_Estimation_with_an_Event_Camera_CVPRW_2024_paper.html?utm_source=chatgpt.com)

Gossard et al. · CVPR Workshops 2024

参考事件光流、球面商标事件提取、旋转周期估计。主要用于开发 光流和周期转速 Baseline，与 s-CMax 比较。

3\. [Event-based Ball Spin Estimation in Sports](https://openaccess.thecvf.com/content/CVPR2024W/CVsports/html/Nakabayashi_Event-based_Ball_Spin_Estimation_in_Sports_CVPRW_2024_paper.html?utm_source=chatgpt.com)

Nakabayashi et al. · CVPR Workshops 2024

研究球面三维投影、事件运动补偿与 CMax 优化。与 Sony s-CMax 的思路最接近，适合用于验证旋转几何、运动补偿和目标函数。

4\. [A Unifying Contrast Maximization Framework for Event Cameras](https://openaccess.thecvf.com/content_cvpr_2018/html/Gallego_A_Unifying_Contrast_CVPR_2018_paper.html?utm_source=chatgpt.com)

Gallego et al. · CVPR 2018

CMax 的基础论文。重点理解 Image of Warped Events（IWE）、运动补偿以及为什么最大化事件图像对比度能够用于运动参数估计。

上述四篇论文分别承担主算法、事件前端、对照算法和理论基础的作用。

[image](https://www.google.com/s2/favicons?domain=https://arxiv.org\&sz=32)

arxiv.org+3



## 二、Sony 官方公开伪代码（非常值得给 Codex）

5\. [Sony Ace — Supplementary Material / Approximate Python Code](https://sonyresearch.github.io/ace_public/pseudo_code/?utm_source=chatgpt.com)

这是 Sony 《Outplaying Elite Table Tennis Players with an Autonomous Robot》的补充材料，包含：

- `events2d_to_3d`：二维事件映射到单位球面。
- `warp_events3d_by_spin`：使用候选角速度补偿事件。
- `event3d_to_spherical`：三维点转换成球面坐标。
- `spherical_to_image_of_warped_events`：生成球面 IWE。
- `SpinRefinerCMax`：通过对比度比较候选旋转参数。

该页面明确说明代码是用于解释关键算法的近似实现，省略了生产系统细节。

Codex 必须注意：官方伪代码本身存在明显笔误，不能直接复制运行。 例如球面深度计算出现 `x3d**2 + x3d**2`，而 Sony 的正式数学公式应为 `x²+y²`；伪代码还存在变量拼写和参数接口不完整等问题。应以正式论文公式为准，并通过单元测试验证。

[image](https://www.google.com/s2/favicons?domain=https://sonyresearch.github.io\&sz=32)

ace_public+2



## 三、DVSense 官方 SDK 与采集文档（P0–P2、P7）

这部分是 Codex 在 Ubuntu 上连接 G1、读取 RAW、建立数据处理工具必须查阅的资料。

| 编号 | 官方资料                                                                                                        | 对应任务                          |
| -- | ----------------------------------------------------------------------------------------------------------- | ----------------------------- |
| 6  | [软件与 SDK 安装](https://sdk.dvsense.com/zh/html/install_zh.html?utm_source=chatgpt.com)                        | P0：Insight、Python SDK、C++ SDK |
| 7  | [相机连接准备](https://sdk.dvsense.com/zh/html/tutorial_electrical_connection_zh.html?utm_source=chatgpt.com)     | P0：G1 千兆以太网连接要求               |
| 8  | [DVSense Insight 上位机操作手册](https://sdk.dvsense.com/zh/html/camera_viewer_manual.html?utm_source=chatgpt.com) | P1、P7：录制、回放、ROI、文件导出          |
| 9  | [RAW 文件读取教程](https://sdk.dvsense.com/zh/html/file_zh.html?utm_source=chatgpt.com)                           | P1：EVT3.0 解码、时间窗读取            |
| 10 | [事件相机及数据格式简介](https://sdk.dvsense.com/zh/html/intro_to_dvs.html?utm_source=chatgpt.com)                     | P1：事件坐标、时间戳、极性                |
| 11 | [相机工具配置](https://sdk.dvsense.com/zh/html/tools_zh.html?utm_source=chatgpt.com)                              | P2、P7：Bias、ROI、滤波和事件速率        |
| 12 | [官方光学参数计算器](https://sdk.dvsense.com/zh/html/cal_opt_param.html?utm_source=chatgpt.com)                      | P7：球的成像直径、工作距离与镜头选择           |

这些页面适用于 DVSLume-G1；其中 RAW 文件读取教程明确说明 G1 录制 `.raw`，采用 EVT3.0 编码。Insight 可以录制和回放，并导出 MP4 预览。

[image](https://www.google.com/s2/favicons?domain=https://sdk.dvsense.com\&sz=32)

sdk.dvsense.com+4



### 13. 官方 Python 示例库

[dvsense_driver_python_sample](https://github.com/dvsense/dvsense_driver_python_sample?utm_source=chatgpt.com)

建议 Codex 重点阅读其中三个文件：

| 源文件                                                                                                                                                                           | 用途                |
| ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------- |
| [raw_file_reader.py](https://github.com/dvsense/dvsense_driver_python_sample/blob/main/raw_file_reader.py?utm_source=chatgpt.com)                                             | RAW 文件加载、读取事件和时间戳 |
| [raw_file_recorder_video_export_sample.py](https://github.com/dvsense/dvsense_driver_python_sample/blob/main/raw_file_recorder_video_export_sample.py?utm_source=chatgpt.com) | RAW 与视频导出         |
| [config_camera_params.py](https://github.com/dvsense/dvsense_driver_python_sample/blob/main/config_camera_params.py?utm_source=chatgpt.com)                                   | 相机参数查询和配置         |

该官方示例仓库标注 MIT 许可证，提供 `RawFileReader` 的调用范例。

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

dvsense/dvsense_driver_python_sample · GitHub



特别提醒 Codex：必须验证 SDK 时间单位。 官方 Python 示例中的 `acc_time` 注释和底层时间戳定义值得核对，不能未经实测就把传入的 `10000` 当成 10 ms 或 10 s。建议用实际返回事件的时间戳差分来确认。

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub+1



## 四、算法开源代码与数学实现（P3–P6）

| 编号 | 项目或技术文档                                                                                                                                                     | 用途                            |
| -- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------- |
| 14 | [TimoStoff/event_utils](https://github.com/TimoStoff/event_utils?utm_source=chatgpt.com)                                                                    | 参考事件运动补偿、IWE 构建、CMax 目标函数及优化器 |
| 15 | [Event-Based Visual Flow — Benosman et al., 2014](https://pubmed.ncbi.nlm.nih.gov/24807038/?utm_source=chatgpt.com)                                         | P4：事件局部时空平面拟合与光流初值            |
| 16 | [SciPy Rotation.from_rotvec](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.from_rotvec.html?utm_source=chatgpt.com) | P4：SO(3) 旋转向量与三维运动补偿          |
| 17 | [SciPy Nelder-Mead](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-neldermead.html?utm_source=chatgpt.com)                                    | P4：无梯度局部优化                    |
| 18 | [SciPy Gaussian Filter](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter.html?utm_source=chatgpt.com)                     | P4：球面 IWE 平滑                  |
| 19 | [SciPy Periodogram](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.periodogram.html?utm_source=chatgpt.com)                              | P3：事件率的候选周期与频率分析              |

`event_utils` 已经提供多种事件运动模型、对比度优化方法和可视化工具，但其 `pure_rotation_warp` 是图像旋转模型，不是 Sony 的完整三维球面 s-CMax；球面补偿仍需独立实现。

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub



光流部分建议优先研究 Benosman 的事件时空平面方法，因为 Sony §3.3.1 明确使用这一类初始化思路。

[image](https://www.google.com/s2/favicons?domain=https://pubmed.ncbi.nlm.nih.gov\&sz=32)

PubMed+1



### 可选：用于合成数据和实验标定

| 编号 | 资料                                                                                                | 用途                                         |
| -- | ------------------------------------------------------------------------------------------------- | ------------------------------------------ |
| 20 | [SensorsINI/v2e](https://github.com/SensorsINI/v2e?utm_source=chatgpt.com)                        | P6：从高时间分辨率图像序列模拟事件相机数据；比手工添加随机事件更接近真实传感器行为 |
| 21 | [AprilTag 3](https://github.com/AprilRobotics/apriltag?utm_source=chatgpt.com)                    | P7：如果使用电机旋转装置，可以辅助确定电机轴相对于参考相机的方向          |
| 22 | [Event-Based Vision: A Survey](https://doi.org/10.1109/TPAMI.2020.3008413?utm_source=chatgpt.com) | 事件相机原理、事件光流、事件表示及视觉算法的系统学习资料               |

v2e 是事件相机模拟器，但不是乒乓球专用旋转数据集。第一版可以先自行构建球面标记旋转的理想事件生成器，等需要增加噪声、像素响应差异时再引入 v2e，避免增加不必要的依赖。

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub+2



## 五、实验设置与质量控制资料（P5、P7）

这里再补充两份与真实数据质量直接相关的官方文档。

23\. [DVSense Bias Configuration](https://sdk.dvsense.com/en/html/bias_en.html?utm_source=chatgpt.com)

G1 使用 IMX636 对应的 Bias 参数组。主要用于理解事件阈值、噪声和事件率之间的关系。官方强调参数范围应以设备实际返回值为准。

[image](https://www.google.com/s2/favicons?domain=https://sdk.dvsense.com\&sz=32)

sdk.dvsense.com



24\. [DVSense 事件信号处理](https://sdk.dvsense.com/zh/html/event_signal_processing_zh.html?utm_source=chatgpt.com)

覆盖防闪烁、事件轨迹过滤以及事件速率控制。对阶段一尤其重要：不要在没有记录配置的情况下启用强过滤，否则可能把有用的周期性事件一并删除。

[image](https://www.google.com/s2/favicons?domain=https://sdk.dvsense.com\&sz=32)

sdk.dvsense.com



此外，Sony 在 §4.2 固定球实验中提供了很好的实验设计参考：

- 使用电机编码器测转速，AprilTag 辅助确定转轴。
- 覆盖不同转速、旋转轴和商标朝向。
- 默认使用 10 ms 事件窗口。
- 对无足够球面特征事件的窗口进行筛除。

论文给出的窗口长度、质量阈值和精度结果都不能直接当作 G1 的最佳参数，需要用自己的数据重新验证。

[image](https://www.google.com/s2/favicons?domain=https://arxiv.org\&sz=32)

arxiv.org



## 六、这些资料如何对应 Codex 的开发任务？

| 任务               | 优先参考编号                      |
| ---------------- | --------------------------- |
| P0：环境准备          | 6、7、13                      |
| P1：RAW 解码与预览     | 8、9、10、13                   |
| P2：球体 ROI、事件筛选   | 1、8、11、23、24                |
| P3：周期转速 Baseline | 2、19                        |
| P4：光流与球面 s-CMax  | 1、3、4、5、14、15、16、17、18      |
| P5：质量评估和报告       | 1（§4.2）、2、3                 |
| P6：合成数据和测试       | 16、17、20                    |
| P7：固定球实验与真值      | 1（§4.2）、8、12、21、23、24       |
| P8：CLI 与验收       | 官方 SDK API 文档，以及项目自身定义的验收标准 |

需要区分： P5 的 JSON 字段、有效窗口比例，以及 P6 的合成测试误差阈值，主要是我们为工程可验证性提出的设计要求，不是 Sony 论文原样提供的接口或标准。

## 七、给 Codex 的补充要求

建议把下面这段追加到之前的开发任务书：

请将我提供的阶段一参考文献及官方链接整理到 `docs/references.md`，按 P0–P8 建立来源映射。

优先精读 Sony ECCV 2026 论文 §3.3、§4.2，以及 Sony Ace 官方公开的 ContrastMaximization Python 伪代码，并与 Gossard 2024、Nakabayashi 2024 进行方法对照。

要求：

1. 每个核心算法模块标明参考论文的章节、公式或公开代码位置。
2. 优先使用 DVSense 官方 RAW Reader，不重新实现已有解码功能。
3. Sony 公开伪代码仅作参考，检查其中的数学、变量和接口错误，以正式论文为准。
4. 论文直接支持的实现、依据论文推导的实现、我们自行选取的工程参数必须明确区分。
5. 记录使用的代码库版本、许可证及其实际复用范围。
6. 不引入本阶段用不到的神经网络、ROS 2 或飞行球追踪模块。

最终交付的 `algorithm.md` 应让审查者能够从代码追溯到论文依据，并能通过数学单元测试验证主要公式。

推荐阅读顺序：Sony ECCV 2026 → Sony 公开伪代码 → Nakabayashi 2024 → Gallego 2018 → Gossard 2024。

这五份材料足够支撑阶段一的核心算法设计。其他资料主要解决相机接入、代码实现和实验验证问题。