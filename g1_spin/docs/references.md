# Phase 1 参考资料与实现溯源

核查日期：2026-10-29。链接为参考资料，不代表全部已集成。算法代码为依据数学描述独立实现；未复制 Sony、event_utils、v2e 或 AprilTag 源代码。依赖实际安装版本应以交付环境的依赖清单为准。

## 核心材料与核查范围

| 材料 | 本次核查及用途 |
| --- | --- |
| [Sony, Event-based Gaze Control System, arXiv 2606.26780v2](https://arxiv.org/html/2606.26780v2) | 已读 §3.3.1、§3.3.2、§4.2；式(2)球面提升、式(3)IWE、式(4)方差；光流投影、大小搜索和 Nelder–Mead；编码器/AprilTag 验证。HTML 未包含其引用的补充材料。 |
| [Sony Ace 公开近似 Python](https://sonyresearch.github.io/ace_public/pseudo_code/) | 已读 `ContrastMaximization`、`SpinRefinerCMax` 和调用者。属于 Nature Ace 补充材料，不是可直接运行的离线 s-CMax 项目；仅交叉核验概念，错误清单见 algorithm.md。 |
| [Gossard 2024](https://arxiv.org/html/2404.09870v1) / [CVF 正式页](https://openaccess.thecvf.com/content/CVPR2024W/CVsports/html/Gossard_Table_Tennis_Ball_Spin_Estimation_with_an_Event_Camera_CVPRW_2024_paper.html) | 已读 arXiv v1 §3.3 式(5)/(6)、事件率周期和 §4.5 局限。周期和光流对照依据；当前周期实现采用工程频谱方案，未逐行复现其 EMA 零交叉法。 |
| [Nakabayashi 2024](https://openaccess.thecvf.com/content/CVPR2024W/CVsports/html/Nakabayashi_Event-based_Ball_Spin_Estimation_in_Sports_CVPRW_2024_paper.html) | 已读取 CVF 全文（页码 3367–3375）；§4.2 式(9)负深度、§4.3 式(10)/(11)四元数逆时补偿、§4.4 丢弃深度投回图像平面；§5.1 使用 50000 事件和 Adam 1500 次迭代。与 Sony 球面 IWE 区分；本项目不声称复现该方法。 |
| [Gallego 2018](https://openaccess.thecvf.com/content_cvpr_2018/html/Gallego_A_Unifying_Contrast_CVPR_2018_paper.html) | CMax/IWE 理论来源，本次具体实现公式以已核验 Sony §3.3.2 为准；未导入其代码。 |
| [Benosman 2014](https://pubmed.ncbi.nlm.nih.gov/24807038/) | Sony 引用的事件平面拟合光流理论来源；已实现工程局部时间平面标量法向约束 Baseline；未声称逐行复现 Benosman 或 Sony 的原始实现。 |

## P0–P8 来源映射

| 阶段 | 资料 | 实际使用与界限 |
| --- | --- | --- |
| P0 | [SDK 安装](https://sdk.dvsense.com/zh/html/install_zh.html)、[连接](https://sdk.dvsense.com/zh/html/tutorial_electrical_connection_zh.html)、官方示例库 | 实际机器及安装接口优先；不修改网络或固件。 |
| P1 | [RAW](https://sdk.dvsense.com/zh/html/file_zh.html)、[Insight](https://sdk.dvsense.com/zh/html/camera_viewer_manual.html)、[事件格式](https://sdk.dvsense.com/zh/html/intro_to_dvs.html)、[reader 示例](https://github.com/dvsense/dvsense_driver_python_sample/blob/main/raw_file_reader.py) | 使用已有官方 SDK reader；已通过一份 840449 事件、1280×720、3.038474 秒 RAW 的检查与 100 ms 预览；未自行解码 EVT3.0，原文件 hash 不变。 |
| P2 | Sony 式(2)、[相机工具](https://sdk.dvsense.com/zh/html/tools_zh.html)、[信号处理](https://sdk.dvsense.com/zh/html/event_signal_processing_zh.html) | 手工圆 ROI；内缩比例、热点阈值是工程参数，禁止默认强滤波破坏商标轨迹。 |
| P3 | Gossard §3.3、[periodogram](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.periodogram.html) | 事件率周期候选；不能证明一周期等于一圈。 |
| P4 | Sony §3.3、Ace `ContrastMaximization`、Nakabayashi、Gallego、Benosman；[Rotation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.from_rotvec.html)、[Nelder–Mead](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-neldermead.html)、[Gaussian](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter.html) | 几何与目标遵循 Sony；多初值粗搜索为工程替代；独立光流 Baseline 已有秩/条件数/残差检查，尚未用于初始化。 |
| P5 | Sony §4.2、Gossard §4.5 | 质量门槛、JSON/CSV 字段、有效比例及一致性检查是项目设计，不是论文置信度定义。 |
| P6 | SciPy Rotation、[v2e](https://github.com/SensorsINI/v2e) | Rotation 独立生成轨迹供验证；v2e 未引入。当前是理想球面点轨迹，不是物理事件模拟。 |
| P7 | Sony §4.2、Insight、[光学计算](https://sdk.dvsense.com/zh/html/cal_opt_param.html)、[AprilTag](https://github.com/AprilRobotics/apriltag)、[Bias](https://sdk.dvsense.com/en/html/bias_en.html) | 参考独立编码器和轴标定设计；AprilTag 未集成，G1 Bias 以设备范围为准。 |
| P8 | 官方 SDK、本任务书的 CLI 和验收要求 | synthetic 与实际 RAW 读取分别验收；已验证 RAW 读取不等于验证真实旋转估计精度。 |

补充原理资料：[Event-Based Vision Survey](https://doi.org/10.1109/TPAMI.2020.3008413)。

## 代码版本、许可与实际复用

以下 HEAD 是核查时远端状态，不表示代码复制自该提交。

| 项目 | 核查版本/许可 | 实际复用 |
| --- | --- | --- |
| [DVSense Python 示例](https://github.com/dvsense/dvsense_driver_python_sample) | HEAD `923943401c8f8d8a5581a6d559a6d7b1e48452ca`；根 LICENSE 已核验为 MIT，copyright 2025 dvsense | 参考 reader 调用；运行使用 Ubuntu 既有 `scripts/event_io.py` 与官方 `dvsense_driver` 1.3.0。示例 MIT 不自动等于二进制 SDK 的许可；SDK 按厂商条款。 |
| [SonyResearch/ace_public](https://github.com/SonyResearch/ace_public) | HEAD `57ac16a42c25c47bde24cd2a592804f4b800085b`；所读页面/README 未见明确代码授权，根常见 LICENSE 路径未取得文件 | 只读算法交叉核验，未复制代码；不将公开可读推定为 MIT。 |
| [TimoStoff/event_utils](https://github.com/TimoStoff/event_utils) | HEAD `dc0a0712156bb0c3659d90b33e211fa58a83a75f`；根 LICENSE 已核验 MIT，copyright 2020 Timo Stoffregen | 阅读 README 确认事件工具和对比度优化范围；未安装、未复制；不能用通用图像旋转代替球面物体旋转。 |
| SciPy / NumPy | 运行依赖；BSD-3-Clause；精确版本见环境清单 | 优化、滤波、数值计算；合成生成器用 SciPy Rotation，估计器 Rodrigues 独立实现可交叉测试。 |
| v2e / AprilTag | 未引入、未锁定版本、未核验本次具体提交许可 | 仅未来模拟/标定选项，不属于当前交付依赖。 |

论文网页的开放访问不等于代码复用许可；CVF 页面明确版权由作者或相应权利人保留。本项目引用论文思想和自行推导公式，没有分发论文正文或第三方源码。
