# G1 固定球三维旋转 Demo

## 2026-11-13 RGB90同步录制与旋转方法对照

本次关闭D455投射器，以640×360 RGB约90FPS与事件相机近似同步录制。
[分析报告](dual_capture/analysis_20261113/report.md)记录文件质量、初步时间关联和RGB周期参考。
[方法汇总](dual_capture/analysis_20261113/methods_comparison.md)比较固定轴、10/50/250ms自由三维与50ms估轴＋250ms估速组合，包含轴向、RPM、rad/s及相对参考偏差。
10ms的1200/3500RPM搜索上限对照结果接近。组合在四个探索窗口的速度偏差约0.3%–1.5%，轴相对粗参考约4–8°，仍需留出片段验证。
录制与分析脚本位于dual_capture；原始录像和大型数组按仓库忽略规则保留在采集机或本地。

最新[斜轴偏高分析](experiments/oblique_bias_20261113/report.md)：发现一维/三维补偿参考时刻不一致；统一后run3一维中位283.0→274.8RPM，三维271.9RPM。两者仍受半径/球心及事件观测影响，原一维汇总保留为历史诊断。

## 2026-11-13 Sony方案与路线裁决

[系统审查与可证伪实验](experiments/strategy_audit_20261113/report.md)确认旧事件含可预测的旋转相容运动信息，
但中心平移与横向自转存在强几何耦合：不同RPM候选可有近似相同的后段像素预测。
保留s-CMAX为条件性拟合器；后续优先固定轴的一维速度、独立中心约束和未来预测验证，
不以继续调网格/轮廓参数替代真实RPM或未知三维轴验收。生产默认未修改。

## 2026-11-13 第一次圆盘录像转速复核

[圆盘转速复核](experiments/disk_recheck_20261113/report.md)完成126条s-CMAX对照和允许平移的纹理复现检查。
222403减速量级约200→130→70RPM在多种中心/窗口假设下仍成立；高平移段补偿使候选改变约10RPM，
但不能证明补偿后绝对转速更准。主数据没有独立旋转真值，保留模型范围而不报告精确真实误差。

## 2026-11-12 新录制前置工具

[固定评测工具入口](validation/pre_capture_20261112/README.md)包含三方案快照、统一离线运行、同步脉冲拟合与参考时间转换。
[前置工作验收](validation/pre_capture_20261112/RESULTS.md)记录8项测试、旧实录回归及84次受控几何误差实验；
[采集清单](validation/pre_capture_20261112/CAPTURE.md)明确同步、跨相机投影和独立录制测试的要求。
这些工具不修改生产默认，不代表已完成真实相机标定或同步精度验收。

## 2026-11-10 电钻近似参考与窗口对照

用户确认60fps手机视频为另拍的同档位参考，可用于当前近似转速对比。
两段贴片过线读数约249.44/244.39 RPM；不作为同步逐窗转速或三维轴真值。
本轮126条离线对照及结论见[球心与窗口实验报告](experiments/drill_window_geometry_20261110/report.md)。
20ms + 三点中值平滑中心暂作候选，生产默认未修改；此前报告中的“无真值”保留其当时的实验口径。
后续[44区间扩展验证](experiments/drill_candidate_validation_20261110/report.md)显示转速偏差与轴波动改善，
但两个基线可用窗因20ms几何缺口丢失，斜轴二维杆投影差增大；候选不推广为默认。
随后[几何缺口与斜轴诊断](experiments/drill_geometry_diagnosis_20261110/report.md)完成48条六组对照；
仅按几何覆盖回退10ms的实验路由在既有44区间回放中恢复两窗，旋转有效24/32。
斜轴仍对球心补偿敏感，统一取消补偿会导致平移误报；此路由未接入生产默认。
进一步[原始轮廓与外杆交叉检查](experiments/drill_center_evidence_20261110/report.md)核验了144个快门：
斜轴旧中心横向变化比杆参考更大，两侧轮廓中心不一致，支持存在拟合干扰，但尚不能唯一分离真实振动与观测误差。
新增[外侧支撑弧候选](experiments/drill_outer_contour_20261110/report.md)：5ms斜轴共同支持下两侧中心差
3.39→1.14px，三个共同旋转窗二维杆投影差中位1.92→0.51度；但丢失一窗且转速未同步改善，仍仅作离线候选。
[外侧前端扩展与回退验证](experiments/drill_outer_validation_20261110/report.md)覆盖18个其他片段：
新旧均8/12旋转可用、纯平移0/4非零；参考偏差13.50→13.33RPM且不同录制收益不一致，不升级默认。

Ubuntu：`pingpong1@172.16.16.225:~/g1_spin/phase1_demo`。
Mac 源码和结果副本：`/Users/vae/Projects/PingPongRobot/g1_spin`。
独立 Python 3.10 虚拟环境；没有覆盖原有 `~/g1_spin/scripts` 或 RAW 文件。

## 立即演示

在 Ubuntu 执行：

```bash
cd ~/g1_spin/phase1_demo
bash scripts/run_demo.sh outputs/my_demo
xdg-open outputs/my_demo/index.html
```

无图形 SSH 终端只运行前两行；把输出目录复制回 Mac 后打开 `index.html`。
已有实测结果：`outputs/demo_final/index.html`，完整数据见 `benchmark.json` 和各场景的 `result.json`、`windows.csv`。
HTML 报告无需联网；MP4 链接可用系统播放器打开。

演示覆盖 500/1000/2000/3000 RPM、反转、10% 背景噪声、空输入，使用 10 ms 窗口、3 个窗口/场景。
每场景输出事件预览、球面补偿对比、优化曲线、RPM 曲线和报告。
输入仅有事件 x/y/t/p，没有把旋转真值或跟踪关系传给估计器。

**这是理想球面点轨迹 demo，不是物理 DVS 仿真或真实旋转精度验收。** 光流模块会对缺少二维边缘支撑的点轨迹返回无效，当前基准中未恢复有效光流；s-CMax 使用独立多轴网格初值。不能宣称完整复现 Sony 的光流初始化流程。

## 重新安装

已有 Ubuntu 环境可直接使用。新环境需 Python 3.10+：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt
.venv/bin/python -m pip install -e . --no-deps --no-build-isolation
.venv/bin/python -m pytest -q
```

`requirements-lock.txt` 固定本次通用依赖；无系统 Python 包继承。官方 SDK 是单独的 Linux Python 3.10 wheel：

```bash
.venv/bin/python -m pip install ~/g1_spin/wheels/dvsense_driver-1.3.0-cp310-cp310-manylinux_2_28_x86_64.whl
```

RAW 路径另外依赖原工程官方 C++ SDK 辅助程序 `~/g1_spin/scripts/sdk_check` 和 `scripts/event_io.py`，解决 1.3.0 的 EOF 尾时间字问题；不会自己解码 EVT3。
其他机器需按原工程 `docs/sdk_validation.md` 安装/编译，或设置 `G1_LEGACY_ROOT` 指向对应工程。没有 SDK 时仍可完整运行合成 demo。
若官方 PyPI 下载超时，可在 pip 命令添加 `-i https://pypi.tuna.tsinghua.edu.cn/simple`。

## 真实 RAW 读取与预览

```bash
.venv/bin/python -m g1_spin.cli inspect ../data/raw/G1_CONNECTIVITY_TEST_20261029_210504.raw --output outputs/inspect
.venv/bin/python -m g1_spin.cli preview ../data/raw/G1_CONNECTIVITY_TEST_20261029_210504.raw --window-ms 2 --duration-ms 100 --output outputs/preview.mp4
```

已实测完整读取 840449 事件，无倒序，生成 50 帧预览并重新解码检查。连通性 RAW 没有已知旋转真值，因此未对它报告球转速。
SDK 可能维护旁路索引文件；RAW 本体 SHA256 在本次读取前后相同。

人工从预览测量 `(cx, cy, radius)`，填入 `configs/phase1.yaml`；三个值默认空，避免误用合成 ROI。
支持配置窗口、步长、绝对 `start_us`、`duration_ms` 和搜索 `max_rpm`：

```bash
.venv/bin/python -m g1_spin.cli analyze ../data/raw/run001.raw --config configs/phase1.yaml --output outputs/run001
```

默认仅分析 100 ms，确认 ROI 和质量后再扩大范围。`duration_ms: null` 读至结束；事件不足、退化或不一致结果以 null 输出。
预览固定 50 fps，2 ms 窗对应 0.1 倍速；时间和帧数写在同名 JSON。
内部 NPZ 是小规模 demo 交换格式，整份载入内存；真实 RAW 按块读，不整份载入。

## 验收与边界

- [实测验收](docs/acceptance.md)：测试命令、结果、待完成项。
- [算法](docs/algorithm.md)：坐标、公式、工程差异、质量门限。
- [参考资料](docs/references.md)：P0–P8 映射、伪代码笔误、版本/许可证/复用范围。
- [录制流程](docs/experiment_protocol.md)：静止负对照、转速标签、重复实验、元数据。

下一步人工操作：固定有非对称标记的球；用独立编码器/转速计记录真值；录制静止与不同转速/轴向各三次；填写 ROI、实验条件和真值来源。不可从本算法输出反填真值。

## 2026-10-29 真实手动圆盘录像（Mac本地）

3段原始RAW和官方SDK解码数据已复制并逐文件校验，完整清单见 `data/archive_manifest.json`。
最新17.58秒录像的分析入口：`outputs/real_222403/index.html`，详细说明见同目录 `report.md`。
另两段的筛选见 `outputs/older_recordings_qc/report.md`。本次候选转速没有独立真值；球心摆动、遮挡和几何误差使三维轴暂不可靠。

Mac 已安装独立 `.venv`，运行无需连接 Ubuntu：

```bash
.venv/bin/python scripts/analyze_recording_222403.py
.venv/bin/python scripts/build_real_report.py
```

`io.source` 现在可直接读取官方SDK导出的 `manifest.json` + 分块NPY，不要求Mac安装官方SDK。
