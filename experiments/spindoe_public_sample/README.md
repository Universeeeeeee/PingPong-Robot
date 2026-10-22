# 乒乓球旋转检测算法与实验

目的：先用官方真实采集的球体图像检查 SpinDOE 能否运行、点阵预测是否与图像一致。该实验不是 G1 事件相机测试，也不提供独立旋转精度认证。

## 算法入口

| 内容 | 实现 | 说明 |
|---|---|---|
| SpinDOE 原版复现及排序/符号诊断 | `run_sample.py` | 通过固定提交下载官方源代码、模型和22帧样本 |
| 多基准点与一一对应匹配 | `point_matcher.py`、`run_improved.py` | [改进记录](IMPROVEMENTS.md) |
| 符号修正与有界圆周相位拟合 | `corrected_spin.py` | 单位为观测坐标系 rad/s，搜索范围不保证排除混叠 |
| 商标暗纹理检测、LK跟踪与球面旋转拟合 | `logo_frontend.py`、`run_logo_video.py` | [前端说明](LOGO_FRONTEND.md)，无已知圆点模板 |
| 已知旋转、丢帧、坐标和失效条件测试 | `test_*.py` | 17项测试；含合成三维文字球，不能替代实拍精度验证 |

仓库保存算法、测试、复现清单和说明。大模型、第三方视频、虚拟环境及本地运行产物不随本次源码上传；上游源代码/模型/样本通过下述脚本按哈希校验下载。运行后产生的 `outputs/` 路径是生成结果。

实现边界：事件光流与球面对比度最大化（s-CMax）目前仅完成方案调研，尚无本地实现；这里的商标前端处理普通图像帧，不是 G1 原始事件解码器。

## 输入与来源

- [官方仓库](https://github.com/cogsys-tuebingen/spindoe)，固定提交 `ebc8a9eac32e42cdbe2978767ef0134d67dfc123`，Apache-2.0。
- `upstream/data/test/`：22 张球体裁剪 PNG，文件名为纳秒时间戳；序列跨度 65.573038 ms。帧间隔中位数约 2.623 ms，存在双倍间隔，不能用固定播放 FPS 替代采集时间戳。
- `upstream/data/model/spindoe.ckpt`：官方预训练点检测器。
- `upstream/cad/ref_pattern.csv`：与样本匹配的已知球面点阵。
- 原文件及 Git blob 哈希记录在 `source_manifest.json`。下载和运行时均核对哈希，未修改上游文件。

官方 `data/trajectories` 是位置与 SpinDOE 角速度结果，不是更多图像，也不是独立旋转真值。

## 运行

进入本目录，使用 Python 3.13 创建独立环境。依赖版本记录于 `requirements.txt`；预览字幕使用 Pillow 自带字体，不依赖 macOS 字体路径。

```sh
python3 download_sample.py
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run_sample.py
.venv/bin/python run_improved.py
.venv/bin/python -m unittest test_logo_frontend test_point_matcher test_corrected_spin -v
```

使用自己的商标视频与球圆标注：

```sh
.venv/bin/python run_logo_video.py input.mp4 --circles circles.json --out outputs/my_logo
```

只有确认采集帧率后才增加 `--capture-fps`；否则仅输出每视频帧角位移，不把播放帧率当测量时钟。`LOGO_FRONTEND.md` 中的本地研究素材回放命令需要另外准备对应输入，并非仓库自带实拍测试集。

此脚本先直接调用上游 `SpinDOE.estimate()`，保留点检测、点阵匹配和 RANSAC 估计。绕过默认 CUDA 构造、GUI 调试入口，并在 CPU 上以 `weights_only=True` 加载官方模型参数。运行副本仅删除未使用的 `dotdataset` 和 `wandb` 导入，避免安装训练数据和日志依赖；原始源码不变，没有训练或上传数据。

若默认算法失败，脚本另运行明示的诊断变体：仅将检测点按 blob 大小降序排列，再进入原有几何匹配和 RANSAC。这不是默认算法的成功复现。诊断变体比较 5 个 RANSAC 随机种子。

## 本次实际结果（2026-10-22）

| 检查 | 结果 |
|---|---|
| 原始算法，默认点顺序 | 7/22 帧通过 DOE；RANSAC 未产生整段角速度 |
| 按点大小降序的诊断变体 | 12/22 帧通过 DOE，10 帧被 RANSAC 保留 |
| 诊断候选转速 | 132.9119 rps；5 个随机种子得到相同值 |
| 诊断候选，原符号前向姿态残差中位数 | 72.45° |
| 诊断候选，反号后前向姿态残差中位数 | 6.39° |
| 独立生成姿态的符号检查 | 输入 `[20,30,-40]` rad/s，上游回归输出 `[-20,-30,40]` rad/s |

点排序敏感性另记录在 `outputs/ordering_diagnostics.json`：反序和距中心排序产生约 97.7–98.8 rps 候选，但前向一致性差；按 x 排序无结果。因此只看到一个转速数字不代表估计可靠。没有针对该片段的独立旋转真值，不能认定 132.9 rps 正确，6.39° 也不是旋转轴真实误差。

## 输出

- `outputs/results.json`：采集时间、姿态成功数、RANSAC 保留帧、角速度、转速、方向诊断及运行时间。
- `outputs/baseline_estimates.npz`：默认算法的逐帧热图和姿态。
- `outputs/estimates.npz`：诊断变体的逐帧时间戳、热图、姿态和保留帧。
- `outputs/diagnostic_60x.mp4`：输入、诊断变体原符号和取反角速度预测的三栏对比；60 倍慢放，按原时间间隔保持上一帧，不插值球体运动。末帧额外显示一个中位采样间隔。
- `outputs/contact_sheet.png` / `comparison_frame.png`：静态对比。

上游 `spin_regressor.py` 的角速度回归含一个带 HACK 注释的负号。脚本不会静默修复它，而是同时报告原输出与取反输出的前向姿态残差，并用已知角速度生成的姿态序列单独检查符号。取反只作为明示的诊断结果，不能替代真实角速度验证。

## 结果解释

- DOE 姿态成功和 RANSAC 保留是两个不同计数；都不能视作真实准确率。
- 姿态残差比较的是同一套图像推导出的姿态与其旋转预测，是内部一致性指标。
- 原始角速度在 DOE 观测坐标中：x 向图像右、y 向图像上、可见球面 z>0。未标定到球台或机器人坐标。
- 单段公开样本没有独立转速或轴向标签，所以不报告真实旋转误差。
- 任意反光点布局、仅商标球、事件累积图，都不满足默认点阵和训练图像的相同条件，需要另行适配。

## 其他公开实拍来源

[openspin examples](https://github.com/JamesQFreeman/openspin/tree/main/examples) 提供 `rally_240fps.mov`（约 32.3 MB，240 fps）。其手绘点阵与本样本不同，不能直接套默认 SpinDOE。这里只列来源，未下载或重新分发该视频。
