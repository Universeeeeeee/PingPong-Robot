# 项目任务：DVSLume-G1 固定乒乓球旋转估计 Phase 1

## 0. 项目背景与总目标

目前已有：
- 一台 DVSLume-G1 事件相机及官方 PoE 供电设备。
- Ubuntu 主机 `pingpong1`，Mac 已通过 SSH 连接。
- DVSense Insight 上位机及 DVSense SDK（请先检查安装状态）。
- 后续将使用 Insight 录制固定位置旋转乒乓球的 `.raw` 原始事件数据。

请在 Ubuntu 上开发一套**离线乒乓球旋转估计系统**，实现：

输入：G1 的 EVT3.0 `.raw` 原始事件数据。

输出：
1. 旋转速度（rpm / rps / rad/s）。
2. 三维角速度向量 `omega = [wx, wy, wz]`。
3. 归一化旋转轴 `axis = omega / ||omega||`。
4. 分析时间窗、估计质量与失败原因。
5. 事件可视化、运动补偿前后对比、转速曲线。
6. CSV/JSON 实验结果及简洁的 Markdown 报告。

仅实现**固定相机、固定球心、离线分析**。不开发飞行球跟踪、ROS 2、神经网络、实时控制或多相机融合。

主要参考文献：

Sony：
https://arxiv.org/html/2606.26780v2

重点阅读：
- Section 3.3 Offline spin estimation
- Section 3.3.1 Optical flow spin estimation initialization
- Section 3.3.2 s-CMax spin estimation refinement
- Section 4.2 Offline spin estimation experiments: ball spinner

Gossard 2024：
https://arxiv.org/abs/2404.09870

DVSense 官方 SDK：
https://sdk.dvsense.com/zh/html/file_zh.html

Python RAW 读取示例：
https://github.com/dvsense/dvsense_driver_python_sample/blob/main/raw_file_reader.py

### 开发原则

- 简单优先，Python + NumPy + SciPy + OpenCV 即可。
- 优先在 CPU 上运行，不强制 CUDA。
- 先实现可解释的 Baseline，再实现 Sony s-CMax。
- 不要写庞大的通用视觉框架。
- 论文未提供的具体实现参数应标注为工程选择。
- 对没有足够事件或无法确定旋转轴的数据，返回无效结果，不要虚构转速。
- 每个模块都必须有明确的输入、输出和测试。
- 不修改 Ubuntu 网络配置、SSH、相机固件或无关工程。
- 先检查现有目录，不覆盖已有采集数据和程序。

---

## P0：检查环境并准备工程

通过现有 SSH 连接 Ubuntu，检查：

1. Ubuntu 版本、CPU、可用内存、磁盘空间。
2. Python 版本与现有虚拟环境。
3. `dvsense_driver`、DVSense Insight、SDK 是否可用。
4. `~/g1_spin/` 中是否已有工程和真实 RAW 数据。

建议工程结构：

```text
~/g1_spin/
├── pyproject.toml
├── README.md
├── configs/
│   └── phase1.yaml
├── data/
│   ├── raw/
│   └── synthetic/
├── src/g1_spin/
│   ├── io.py
│   ├── events.py
│   ├── ball_roi.py
│   ├── baseline.py
│   ├── optical_flow.py
│   ├── spherical_cmax.py
│   ├── evaluation.py
│   └── cli.py
├── tests/
├── outputs/
└── docs/
    ├── algorithm.md
    └── experiment_protocol.md
```

模块可以合并，但不要过度抽象。

创建独立虚拟环境。优先使用 NumPy、SciPy、OpenCV、Matplotlib、PyYAML、pytest 和官方 SDK。

**完成标准：** 依赖可安装、包可导入、测试命令能启动，README 能指导他人复现环境。

---

## P1：实现 G1 RAW 读取与事件可视化

优先使用官方 `dvsense_driver.raw_file_reader.RawFileReader`，参考：

```python
from dvsense_driver.raw_file_reader import RawFileReader

reader = RawFileReader(raw_path)
reader.load_file()
events = reader.get_n_time_events(interval)
```

注意：上面只是入口示例，必须检查本机实际 SDK 版本、返回的事件字段、时间戳单位以及时间区间参数单位。

不要自行猜测 EVT3.0 的二进制编码，也不要把官方 RAW 当成普通视频读取。

### 1. 统一事件格式

转换为内部结构：

```python
EventBatch:
    x
    y
    t_us
    polarity
```

要求：
- 保留原始整数微秒时间戳。
- 使用相对时间进行浮点计算，避免大时间戳精度损失。
- 检查时间戳顺序、重复、跳变和异常值。
- 检查坐标是否超出传感器分辨率。
- 事件按时间块读取，避免整份文件占满内存。
- 原始 RAW 只读，不修改。

### 2. 数据检查命令

实现：

```bash
python -m g1_spin.cli inspect data/raw/test.raw
```

输出：
- 设备数据分辨率
- 起止时间戳和录制时长
- 事件总数和平均事件率
- ON/OFF 事件比例
- 时间戳质量检查结果
- 每毫秒事件数量曲线

### 3. 可视化命令

实现：

```bash
python -m g1_spin.cli preview data/raw/test.raw \
  --window-ms 2 \
  --output outputs/preview.mp4
```

功能：
- 支持 1、2、5、10 ms 等事件累计时间窗。
- ON/OFF 事件使用不同颜色。
- 支持导出 PNG 和 MP4。
- 输出事件计数热力图和时间表面图。
- 允许设置时间范围，跳过没有意义的空白片段。
- 无 GUI 环境也能够导出文件。

**完成标准：** 真实 RAW 文件能解析并回放；在没有真实 RAW 时，使用模拟事件检查完整处理链路。

---

## P2：球体 ROI 与有效事件提取

本阶段球心基本不移动，不需要开发复杂的球检测与追踪系统。

### 1. 首先实现人工 ROI

允许通过配置指定：

```yaml
ball:
  center_x: null
  center_y: null
  radius_px: null
  inner_radius_ratio: 0.9
```

提供预览图，并支持从图像上选择球心和半径；如果 SSH 无图形界面，则允许直接输入三个数值。

注意：

固定球自转时，外轮廓可能没有足够事件，不能假设霍夫圆每次都能自动检测成功。

### 2. 事件筛选

保留球体内部事件：

- 丢弃球体外部背景事件。
- 默认忽略靠近球轮廓的部分事件，降低边缘几何误差。
- 检查热点像素和异常噪声。
- 不使用会严重破坏商标运动轨迹的强滤波。
- 不把所有边缘事件误认为商标事件。

保存：
- 原始事件累积图
- 球体 ROI 示意图
- ROI 内外事件数量统计
- 有效事件率随时间变化曲线

**完成标准：** 能清楚显示被用于旋转估计的事件分布，球心和半径可以重复使用。

---

## P3：实现简单转速 Baseline

不要一开始直接实现复杂的三维优化。

先做一个用于对照的频率估计方法：

1. 统计球体内部的事件率随时间变化。
2. 使用自相关或频谱方法寻找周期性。
3. 估计候选旋转周期 `T`。
4. 转换为候选转速：

```text
rps = 1 / T
rpm = 60 / T
```

注意：
- 事件周期不一定等于球旋转一圈的周期。
- 一个球上存在多个相似标记时，可能发生倍频。
- 单个商标在某些轴向上可能始终可见或完全不可见。
- 背景闪烁或设备噪声也可能出现周期性。

因此该方法只能输出**带条件的候选转速**，不能保证就是球的真实转速。

必须报告周期强度、候选频率和可能的倍频歧义。

**完成标准：** 能对合成周期信号给出正确结果，并在无周期、噪声或多重周期情况下报告不确定性。

---

## P4：复现 Sony 的离线三维旋转估计

这是项目核心。

### 4.1 坐标与几何定义

从球体 ROI 中取得事件位置：

```text
e = (u, v, timestamp, polarity)
```

按照 Sony Section 3.3.2 的正交投影近似，把图像事件映射到单位球面：

```text
x = (u - cx) / r
y = (v - cy) / r
z = -sqrt(1 - x² - y²)
```

其中，`cx`、`cy` 是球心，`r` 是图像球半径。

要求：

- 固定并文档化相机坐标系和旋转方向约定。
- 确保 `x² + y²` 在合法范围。
- 注意 z 的符号与论文一致。
- 对投影边缘的数值不稳定进行处理。
- 编写单位球面映射测试。
- 不把像素平面旋转直接等同于球体三维自转。

### 4.2 事件光流初值

参考 Sony Section 3.3.1：

- 基于事件时空局部平面拟合估计法向光流。
- 剔除拟合不稳定或样本数量不足的区域。
- 把二维运动约束转换为球面切向运动约束。
- 利用 `v = omega × x` 的关系估计三维角速度初值。
- 通过最小二乘求解粗略 `omega_init`。
- 检查线性系统的秩、条件数和估计稳定性。

必须处理 aperture problem，即局部边缘光流只能提供有限方向约束。

如果真实数据暂时不足以稳定得到光流初值，可以先使用可控的多初值粗搜索实现 s-CMax，然后再补充光流初始化，但必须如实标明与 Sony 原方法的差别。

### 4.3 实现球面事件补偿

对每个候选角速度：

```text
omega = [wx, wy, wz]
dt = (event_time - reference_time)
```

通过 SO(3) 旋转矩阵将事件反向补偿到参考时间：

```text
x_warped = R(-omega * dt) @ x
```

之后将三维球面点转换为球面坐标：

```text
theta = asin(z_warped)
phi = atan2(y_warped, x_warped)
```

对球面坐标进行双线性插值累积，构建 Image of Warped Events（IWE）。

必须注意：
- 时间使用秒，角速度使用 rad/s。
- 旋转矩阵必须正确处理角度和轴。
- `phi` 的周期边界必须正确。
- IWE 是球面坐标累积，不是简单的二维图像旋转。
- 相同事件数量、不同角速度条件下，目标函数应具有可比较性。

### 4.4 实现 s-CMax 优化

目标函数：

```text
omega_best = argmax Var(IWE(omega))
```

按论文思路：

1. 使用光流初值或候选旋转轴。
2. 在可配置的转速范围中搜索旋转大小。
3. 使用 SciPy Nelder-Mead 进一步优化角速度三个分量。
4. 使用双线性事件累积及高斯平滑，参考论文的 5×5 核与 sigma=1。
5. 记录目标函数初始值、最终值、迭代次数和收敛状态。
6. 对多个不同初值检查局部最优问题。

搜索范围必须可配置，不能硬编码为某个固定转速。

绘制：
- 原始球面事件累积图
- 经过运动补偿的球面累积图
- 优化过程的目标函数变化
- 候选转速与目标函数的关系曲线

### 4.5 输出角速度

输出：

```json
{
  "valid": true,
  "omega_rad_s": [0, 0, 0],
  "axis_camera": [0, 0, 1],
  "spin_rad_s": 0,
  "spin_rps": 0,
  "spin_rpm": 0,
  "method": "s_cmax",
  "quality": {},
  "failure_reason": null
}
```

上述数值仅是格式示意，不能作为实际估计结果。

公式：

```text
spin_rad_s = norm(omega)
spin_rps = spin_rad_s / (2*pi)
spin_rpm = spin_rps * 60
axis = omega / norm(omega)
```

对于无旋转证据或无法估计轴的情况，必须使用 `valid=false` 和空值，不要错误地返回 `[0,0,1]`。

---

## P5：质量评估与结果报告

对录制数据使用可配置的短时间窗口分析。

默认先测试 10 ms，并支持 5、20 ms 及可配置滑动步长。不要假设同一个窗口长度适用于所有转速。

每个窗口记录：

- 起止时间
- 有效事件数量
- 球面特征是否可见
- 初始角速度
- 优化后的角速度
- 转速 RPM
- 旋转轴
- 优化目标函数与改进幅度
- 是否收敛
- 估计有效性和失败原因

实现基础质量检查：

1. 有效事件不足。
2. 商标信息不可观测。
3. 光流线性系统退化。
4. 优化不收敛或结果触及参数边界。
5. 多初值产生不一致的结果。
6. 相邻时间窗估计严重不一致。
7. ROI 中存在明显的背景干扰。

需要注意：优化收敛、目标函数较高，并不自动证明结果正确。

### 结果汇总

生成：

```text
outputs/<experiment_id>/
├── result.json
├── windows.csv
├── event_preview.mp4
├── rpm_curve.png
├── spherical_before_after.png
├── optimization.png
└── report.md
```

报告应包含：
- 有效窗口比例。
- 有效窗口 RPM 的中位数及离散程度。
- 旋转轴的估计分布。
- 主要失败原因。
- 数据质量说明。
- 如果有真实转速标签，输出相对误差及绝对误差。

---

## P6：合成数据与单元测试

在真实 G1 RAW 文件未提供之前，先完成以下测试。

### 1. 合成旋转球

创建一个具有非对称表面图案的单位球，设置已知旋转轴、角速度和时间长度，模拟图案运动产生的事件。

至少覆盖：

- 正向与反向旋转。
- 不同的三维旋转轴。
- 500、1000、2000、3000 rpm。
- 不同图案可见程度。
- 添加一定数量的背景噪声事件。
- 没有有效事件的窗口。

模拟事件不能直接把真值角速度或跟踪对应关系泄漏给估计器。

### 2. 数学单元测试

至少覆盖：
- 微秒到秒的时间转换。
- rad/s、rps、rpm 换算。
- SO(3) 正向和逆向旋转。
- 单位球面提升。
- 球面角坐标转换与边界。
- 球面事件双线性累积。
- 无事件、时间戳异常、非法 ROI。
- 已知转速下的估计恢复。

### 3. 算法对照实验

对同一组有效合成数据，分别运行：

A. 周期估计 Baseline

B. 光流估计 Baseline

C. s-CMax 完整估计

比较各方法的：
- RPM 绝对误差。
- 旋转轴角度误差。
- 有效窗口比例。
- 计算时间。

在理想合成条件下，可将 RPM 相对误差小于 5%、旋转轴误差小于 5°作为初步调试目标，不作为未经真实实验验证的性能声明。

对于图案不可观测或严重退化的样本，预期应返回失败，而不是为了提高测试通过率强行给出估计。

---

## P7：准备真实 G1 录制实验

在 `docs/experiment_protocol.md` 中写出具体录制步骤，包括：

### 实验环境

- G1 固定在稳定支架上。
- 乒乓球固定球心并绕轴自转。
- 优先使用具有独立转速参考的旋转装置。
- 球面保留明显、非对称的商标或人工标记。
- 选择稳定且充足的照明。
- 确保球面纹理能产生足够事件。
- 尽量降低背景闪烁、振动与运动干扰。

### 采集设计

建议初步录制：
- 静止球负对照。
- 已知低转速。
- 已知中转速。
- 已知高转速。
- 改变旋转轴方向。
- 分别测试有标记和无明显标记的球。

每种条件至少准备三次重复实验，保存完整原始 RAW。

记录：
- 实验编号。
- 相机型号和固件/SDK 版本。
- 镜头和工作距离。
- 图像中的球半径。
- 光照和 Bias 设置。
- 球面标记类型。
- 参考转速和参考旋转轴（如果有）。
- 原始文件路径、SHA256 和采集时间。

所有转速真值必须注明测量来源，不能根据算法输出反过来填写真值。

---

## P8：CLI 与最终验收

提供统一命令行入口，至少实现：

```bash
python -m g1_spin.cli inspect data/raw/test.raw

python -m g1_spin.cli preview data/raw/test.raw \
  --window-ms 2

python -m g1_spin.cli analyze data/raw/test.raw \
  --config configs/phase1.yaml \
  --output outputs/test

python -m g1_spin.cli synthetic-benchmark \
  --output outputs/synthetic
```

真实输入文件的路径、ROI、时间范围和参数均应可配置。

如果当前没有 RAW 文件：
- 完成依赖安装、合成数据生成、算法实现及单元测试。
- 实际运行 synthetic-benchmark。
- 标注 RAW 输入及真实相机测试为待验证。
- 不宣称完成真实设备测试或达到 Sony 论文精度。

### 最终交付要求

1. 完整可运行的源代码。
2. 可重复安装的环境与依赖文件。
3. 通过的单元测试及实际执行命令。
4. 合成数据分析结果及图表。
5. 简洁的录制操作说明。
6. 清晰的算法实现说明和论文对应关系。
7. 当前问题、未完成项及风险清单。
8. 下一步只需要我执行的人工操作。

请按 P0 → P8 的顺序实施，每个阶段先验证再继续。优先保证 RAW 读取、球面几何和数学正确性，不要为了凑齐功能而牺牲可验证性。

如果遇到 SDK API 不兼容，先阅读本机实际接口和官方示例再修改；如果论文无法提供完整实现细节，请采用最小合理工程方案并明确记录差异，不得声称逐行复现 Sony 未开源代码。

请直接操作 Ubuntu 环境，完成能在没有真实 RAW 数据时完成的所有工作，最终给出可执行的验收结果。