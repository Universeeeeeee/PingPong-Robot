# DVSLume-G1 乒乓球旋转录制准备

工程位于 Ubuntu `pingpong1@172.16.16.225:/home/pingpong1/g1_spin/`。
Ubuntu 22.04.5 / x86_64；Insight、C++ SDK、预编译例程、Python SDK 均为官方 **1.3.0**。
文件中的 `SYNTHETIC` 数据只用于离线验证，不是真实相机采集。

## 环境与启动

```bash
cd ~/g1_spin
source .venv/bin/activate
python -m pip check
```

独立 `.venv` 使用系统 Python 3.10；不依赖、不修改现有 conda base。
依赖固定在 `requirements.txt`，实装版本在 `requirements-lock.txt`；官方 PyPI wheel 已保存在
`wheels/`，校验值见 `wheel_sha256.txt`。需要重建时：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --no-index --find-links wheels -r requirements.txt
```

系统包安装脚本：`sudo bash scripts/install_system.sh`。它添加官方带签名 apt 源并安装指定包；
例程包正确名称是复数 **dvsense-driver-samples**。不执行全系统升级或旧包清理。
安装过程中 apt 为 `python3.10-venv` 配套升级了 7 个 Python 3.10 运行库/开发包补丁；
现有 conda 环境未更改。详情在 `logs/install_system.log`。

**在 Ubuntu 本地图形桌面的终端启动：**

```bash
DvsenseInsight
```

桌面应用菜单名称为 DvsenseInsight，实际二进制为 `/opt/DvsenseInsight/DvsenseInsight`，
`/usr/bin/DvsenseInsight` 是其链接。SSH 本身没有图形转发；已经检查 Ubuntu X11 桌面存在，
程序分别在 offscreen 和该桌面运行至超时，日志打印 1.3.0，未发生启动崩溃。
尚未人工确认窗口操作。首次启动日志有 `Create raw file path failed`，因此录制前必须在
General Setting → File 将保存目录设为 `/home/pingpong1/g1_spin/data/raw`，验证可写。

## G1 接线与采集验收

官方说明 G1 使用**千兆以太网**，不走 USB。已识别可选 RJ45 网口：

| 网卡 | 支持速率 | 当前链路 |
| --- | --- | --- |
| eno4 | 10/100/1000/2500 Mb/s，igc | 无链路 |
| eno6、eno7 | 100/1000/2500/5000/10000 Mb/s，ixgbe | 无链路 |
| eno1np0、eno2np1 | E810 SFP 网口，ice | 无链路；需适配模块，不作为首选 |

SSH 使用 Wi‑Fi **wlo5，172.16.16.225/23，默认网关 172.16.16.1**。
本工程没有更改任何网卡、IP、DNS、默认路由、MTU、防火墙或相机 IP。

1. 按相机随附要求供电，用合格网线接 G1 与 `eno4`（或支持 1GbE 的 eno6/eno7）。
2. 运行 `ethtool eno4`，确认实际 `Link detected: yes`、`Speed: 1000Mb/s`、全双工。
   支持 1GbE 不等于当前已协商 1GbE。必要时对实际所接网卡执行同样检查。
3. 只读发现（不要与 Insight 同时占用相机）：

```bash
cd ~/g1_spin
timeout 15s .venv/bin/python scripts/probe_camera.py --out logs/camera_probe.json
timeout 15s scripts/sdk_check
```

4. 以实际发现的序列号测试 3 秒事件及时间戳，并导出真实相机配置：

```bash
timeout 15s .venv/bin/python scripts/probe_camera.py --serial ACTUAL_SERIAL \
  --seconds 3 --export-config configs/camera_actual.json --out logs/camera_capture.json
```

当前 Python/C++ 发现均为 0 台，网络收包、硬件事件回调、时间戳和实际协商速率均**待连接硬件验证**。
若物理链路正常但找不到设备，先核对设备标牌/厂商工具给出的相机 IP 和子网，再确定独立相机网口配置。
不要套用猜测的默认相机 IP；不要运行会改 IP 的 `dvsense_eth_camera_ip_config`，
也不要触碰 wlo5。网络配置变更需明确实际目标后另行处理。

## 第一次 RAW 录制：需要你在 Ubuntu 桌面操作

1. 打开 `DvsenseInsight`，在 General Setting → File 设置上述 `data/raw` 保存目录。
2. Device → OpenCamera 选择实际 DVSLume-G1；确认实时事件及时间戳，记录实际分辨率与序列号。
3. 复制 `configs/experiment_template.json` 到该次实验的同名 JSON，填写实际设置。
   必填：实验编号、相机/固件、分辨率、ROI、所有 Bias/滤波开关和数值、镜头/光圈、光照、
   球标记、距离。未知值保留 null，不能冒充默认值。也保存相机导出的配置和标记照片。
4. 使用高对比度、不对称的球体标记及稳定光照；先测静止场景噪声，再用 1/2/5ms 显示累积检查运动。
   显示窗仅影响预览；ROI、Bias 和片上滤波/限速可能改变实际保存事件，须记录当前值。
5. 点击录制按钮开始，先保留约 1 秒准备动作，再做短实验；再次点击停止并等待文件关闭。
   正式数据必须保留原始 `.raw`，不要把 MP4 当原始数据。Insight 自动按日期时间命名。
6. File → Open Recording（某些界面拼写为 Open Recodring）打开刚生成的 RAW，确认可回放。
   首次实际文件还必须经过以下统计和分析验收。SDK 的人工 RAW 起始块边界见
   `docs/sdk_validation.md`；准备动作不能代替真实文件完整性验证。

```bash
sha256sum data/raw/ACTUAL_RECORDING.raw
python scripts/g1.py stats data/raw/ACTUAL_RECORDING.raw --out logs/run001_stats.json
python scripts/g1.py decode data/raw/ACTUAL_RECORDING.raw --out data/preview/run001_events
```

解码目录保存分块 `.npy` 和 manifest：`t` 为未重置的 uint64 微秒时间戳，`x/y` 为原坐标，
`p=0/1` 为原极性。目录必须不存在，以免覆盖数据。SDK 自身会先建立文件索引；
本工程之后按批次导出，不将完整 RAW 一次读入 Python 内存。原始 RAW 不改写。
统计检查倒退、相同时间戳、坐标/极性异常、事件速率和大间隔。
相同时间戳是合法现象；大间隔不单独证明丢包。倒退会阻止预览/分析，不自动排序掩盖问题。

## 预览与二维 Baseline

```bash
# 三种事件窗；默认只处理开始后的1000ms，避免意外导出超长视频
for win in 1 2 5; do
  python scripts/g1.py preview data/preview/run001_events \
    --window-ms "$win" --fps 100 --out "data/preview/run001_${win}ms.mp4"
done

# 自动圆检测；建议先从导出的 PNG 测量并检查球体区域
python scripts/g1.py analyze data/preview/run001_events \
  --window-ms 2 --out data/preview/run001_spin.csv

# 已测得静止球中心/半径时可显式指定；以下数字是合成测试例子，真实数据须替换
python scripts/g1.py analyze data/preview/SYNTHETIC_clockwise.npz \
  --window-ms 2 --circle 128 128 70 --out data/preview/SYNTHETIC_spin.csv
```

也可直接以 `.raw` 为输入。`--start-us` 接受原始事件时钟的绝对微秒；
`--duration-ms 0` 处理全部剩余数据。MP4 默认 100fps：1/2/5ms 窗分别以
0.1/0.2/0.5 倍速度回放，速度和事件起止时间写入同名 JSON；每个窗口只对应一个视频帧。
同名 PNG 是首个累积窗。分析时间间隔始终来自事件时钟，与 MP4 播放速度无关。

CSV 包含每窗状态、跟踪内点数、有符号 rad/s 和 rpm；JSON 包含片段中位数和 IQR。
从相机看，图像 x 向右、y 向下，**正号为顺时针**。缺少纹理/可靠跟踪时输出失败状态，
有效对不足时不提供片段角速度。只有二维表观旋转，不能据此声称完整三维自旋；
快速运动可能混叠，球面投影也会导致偏差。三维轴和 OptiTrack 真值对齐留待下一阶段。

最小扩展接口在 `scripts/spin_baseline.py`：`detect_ball`、`track_features`、
`estimate_rotation`、`analyze_frames`。C++ 例子在 `scripts/sdk_check.cpp`。

## 无相机复现

```bash
cd ~/g1_spin
source .venv/bin/activate
g++ -std=c++20 scripts/sdk_check.cpp -o scripts/sdk_check \
  -lDvsenseDriver -lDvsenseBase -lDvsenseHal
python -m unittest discover -s tests -v
python scripts/g1.py synthetic --omega 20 --out data/preview/SYNTHETIC_clockwise.npz
python scripts/g1.py preview data/preview/SYNTHETIC_clockwise.npz \
  --window-ms 2 --out data/preview/SYNTHETIC_clockwise_2ms.mp4
```

日志：`logs/tests.log`、`logs/environment_check.txt`、`logs/camera_probe.json`、
`logs/cpp_discovery.log`、`logs/insight_desktop.log`。
完整现场状态与测试结果见 `docs/verification_report.md`。
