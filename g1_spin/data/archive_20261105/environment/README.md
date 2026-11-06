> 保存目录已修复并完成真实 GUI 录制落盘验证，见 [保存路径修复记录](docs/storage_fix_20261029.md)。桌面“G1录制文件”可直接进入 RAW 目录。


> 最新确认（2026-10-29）：用户在修正 eno4 为 192.168.33.10/24 后，已在 Insight 中看到实时事件画面。静态配置已持久保存并启用自动连接；网关/DNS 留空，默认路由仍走 wlo5。此确认仅覆盖实时显示，不等同于完成本地 GUI RAW 录制/回放验收。证据：logs/insight_gui_verified.json。
# DVSLume-G1 旋转实验

Ubuntu 工程：`pingpong1@172.16.16.225:/home/pingpong1/g1_spin/`。
Insight / C++ SDK / 例程 / Python SDK 均为官方 **1.3.0**；独立 Python 3.10 `.venv`。
**2026-10-29：真实 G1 已发现，C++ 实时取流和录制通过；预热后 RAW 的 Python 离线统计、完整分块解码及 2ms 预览验证通过。**
Python 实时探针仍超时，Insight GUI 实时画面、回放与正式旋转实验仍待验收。
`eno4` 专用静态 IPv4 已持久修复为 `192.168.33.10/24` 并激活；SDK 已重新发现 G1。
Insight 画面正等待用户点击 Device → OpenCamera 后反馈；详见 [本次硬件检查](docs/connected_check_20261029.md)。

## 启动和录制

在 Ubuntu 本地图形桌面的终端运行：

```bash
cd ~/g1_spin
source .venv/bin/activate
DvsenseInsight
```

1. G1 按随附要求供电，接千兆以太网。优先 `eno4`，也可用 `eno6/eno7`。
   `ethtool eno4` 本次显示 `Link detected: yes`、`Speed: 1000Mb/s`、全双工。
   相机为 `192.168.33.34/24`；`eno4` 已配置并激活同网段静态地址 `192.168.33.10/24`。
2. Insight → General Setting → File 设置保存目录 `/home/pingpong1/g1_spin/data/raw`。
   首次启动曾报告默认保存路径创建失败，录制前务必设置并确认可写。
3. Device → OpenCamera 选 G1，确认实时事件；复制 `configs/experiment_template.json`，
   填写实际序列号、分辨率、ROI、Bias/滤波、镜头、光照、标记、距离与实验编号。
4. 开启取流后先预热；本轮等待 10 秒后验证通过，但不保证每次重启等待 10 秒均足够。
   录制保留约 1 秒准备动作及短实验，再次点击停止并保存 `.raw` 原件；
   正式录制须检查无时间戳倒退，准备段与预热不能替代逐文件检查。
   File → Open Recording 打开并回放，随后执行统计检查。

SSH 使用 Wi-Fi **wlo5**，其 IP、DNS 和默认路由未改。初次检测后的临时地址曾撤回；
后续已将 `eno4` 的错误 `/0` 前缀修正为 `/24` 并持久激活，设置自动连接且不提供默认路由，
相机接口网关与 DNS 留空。旧重复连接仅禁用自动连接，未删除。
完整接线、发现、配置导出和录制说明见 [录制指南](docs/recording_guide.md)。

## 离线处理

将下例 `run001.raw` 替换成实际录制文件名：

```bash
cd ~/g1_spin
source .venv/bin/activate
python scripts/g1.py stats data/raw/run001.raw --out logs/run001_stats.json
python scripts/g1.py decode data/raw/run001.raw --out data/preview/run001_events
for win in 1 2 5; do
  python scripts/g1.py preview data/preview/run001_events --window-ms "$win" \
    --out "data/preview/run001_${win}ms.mp4"
done
python scripts/g1.py analyze data/preview/run001_events --window-ms 2 \
  --out data/preview/run001_spin.csv
```

RAW 由官方 SDK 读取；分块 NPY 保留未重置的 uint64 微秒时间戳、坐标和 0/1 极性。
倒退时间戳会阻止分析，不自动排序。程序默认预览/分析前 1000ms；`--start-us` 为绝对
事件时间戳，`--duration-ms 0` 处理全部剩余数据。视频默认100fps，1/2/5ms窗分别以
0.1/0.2/0.5倍速度回放，同名 JSON 记录时间和速度，PNG 保存首窗。

分析支持自动球体圆检测，或以 `--circle CX CY R` 输入实际测量区域；输出 CSV/JSON 中的
rad/s、rpm、跟踪质量和失败状态。正号为相机图像中的顺时针（x向右、y向下）。
**结果仅为二维表观角速度**，三维旋转轴和 OptiTrack 真值对齐属于下一阶段。

## 复现测试与交付文档

```bash
cd ~/g1_spin
bash scripts/run_offline_checks.sh
```

当前 14 项测试及端到端流水线通过：包含 SDK 对合成 RAW 的逐事件比较、C++ 复核、
时间戳检查、1/2/5ms视频重新解码、双向与静止旋转基线。所有测试数据均标为 `SYNTHETIC`。
最近合成流水线输出位置见 `logs/latest_synthetic_validation_path.txt`；
新增末尾纯时间字 EOF 回归后，14 项测试通过记录见 `logs/eof_regression.log`。

- [初装验证报告](docs/verification_report.md)：保留初装时的版本和合成测试事实。
- [2026-10-29 硬件检查](docs/connected_check_20261029.md)：相机发现、真实 RAW 统计及剩余限制。
- [SDK 验证及限制](docs/sdk_validation.md)：人工 RAW 索引边界、真实 RAW 的 EOF 修复与 Python 实时限制。
- [录制指南](docs/recording_guide.md)：完整命令、环境重建和相机检查。

目录：`data/raw` 原始录制，`data/preview` 派生数据，`scripts` 程序，`configs` 元数据，
`logs` 实测日志，`docs` 文档，`tests` 合成测试。依赖版本见 `requirements-lock.txt`；
官方 wheel 保存在 `wheels/`，可离线重建环境，不修改已有 conda。
