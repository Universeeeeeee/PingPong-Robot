# G1 固定球三维旋转 Demo

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
