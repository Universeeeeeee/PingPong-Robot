
> 最新确认（2026-10-29）：用户在修正 eno4 为 192.168.33.10/24 后，已在 Insight 中看到实时事件画面。静态配置已持久保存并启用自动连接；网关/DNS 留空，默认路由仍走 wlo5。此确认仅覆盖实时显示，不等同于完成本地 GUI RAW 录制/回放验收。证据：logs/insight_gui_verified.json。
# G1 实机连接检查（2026-10-29）

工程：`/home/pingpong1/g1_spin`。官方 Insight/C++/Python SDK 版本均为 1.3.0。
本次确认真实相机发现、C++ 打开/取流/录制、预热 RAW 的 Python 离线统计、
完整分块解码及 2ms 预览验证；
未完成 Insight GUI 实时画面、回放、Python 实时取流或真实球体旋转精度验收。
当前 `eno4` 静态网络已持久修复并激活；GUI 等待用户点击 Device → OpenCamera 后反馈。

## 链路与配置

- 接口是 `eno4`（不是 `en04`）：1000Mb/s、全双工、Link yes。
- 初始无 IPv4，NetworkManager DHCP 反复失败。无地址二层 GVCP 只读发现
  结果见 `logs/gvcp_l2_identified.json`。
- 相机 DVSLume-G1，序列号 DSP03F100025，固件 V1.0.0；
  IPv4 `192.168.33.34/24`，网关 `192.168.33.1`，分辨率 1280×720。
- 初次检测期间仅为 `eno4` 临时添加 `192.168.33.10/24` 并短暂暂停 DHCP；
  当时检测结束已删除临时地址并恢复自动连接。这是历史状态，当前静态配置见下节。
- SSH 的 `wlo5` 保持 `172.16.16.225/23`，默认网关 `172.16.16.1`，DNS 未改。
- 实际配置导出为 `configs/camera_actual.json`：ROI disabled、全画幅 1280×720，
  Bias 全 0，AntiFlicker/ERC/Trail 均 false，TriggerIn true。

## 后续静态网络修复（当前状态）

用户反馈已配置 IP 但 Insight 无画面后，发现两个 `eno4` 连接均误填
`192.168.33.10/0`，设备未正确安装 IPv4。已修复用户新建的连接 `eno4`
（UUID `93105e1f-2f5b-469d-a23c-716f71445b30`）：地址
`192.168.33.10/24`、`ipv4.never-default yes`、`connection.autoconnect yes`，
网关与 DNS 留空。旧重复连接“有线连接3”
（UUID `b16c3ec6-8553-3e91-9415-9837885398d1`）仅关闭自动连接，未删除。

已激活并核验设备实际 IPv4 为 `192.168.33.10/24`，具有
`192.168.33.0/24` 直连路由；默认路由仍经 `wlo5` 的 `172.16.16.1`。
SDK 再次发现 1 台 DVSLume-G1，DSP03F100025，`192.168.33.34`。
修复前备份为 `logs/eno4_profile_before_fix.txt` 和
`logs/eno4_duplicate_before_fix.txt`。这些结果确认网络配置和 SDK 发现，
不代表 Insight GUI 画面已验证；等待用户点击 Device → OpenCamera 后反馈。

## 两次真实录制

| 检查 | 启动后立即录 3 秒 | 启动后等待 10 秒，再录 3 秒 |
| --- | --- | --- |
| RAW（位于 `data/raw/`） | `G1_CONNECTIVITY_TEST_20261029_205800.raw` | `G1_CONNECTIVITY_TEST_20261029_210504.raw` |
| C++ 直播事件 | 897419 | 860425 |
| 直播时间戳倒退 | 743 | 0 |
| 直播非法事件 | 0 | 0 |

第一次 RAW 的独立时间字检查确认 TIME_HIGH 从 1824 回到 0，即约 7.47 秒
回零；不是仅凭上层分析判断的异常。原文件保留，不排序、不修补时间戳。

第二次的 `logs/hardware_warmed_raw_stats.json`：840449 个 CD 事件，ON 250114、
OFF 590335；时间 9959239–12997713 us，跨度 3038474 us；时间戳倒退 0，
非法事件 0，最大正向间隔 49 us。总数与官方 C++ `getMaxEvents()` 一致。
SDK 文件结束时间为 12997718 us，尾部时间字使 Python EOF 状态未正常结束；
已改用官方索引总数停止读取，详见 [SDK 限制与修复](sdk_validation.md)。
直播计数与录制统计覆盖区间不同；这些检查不能证明网络无丢包。

## 离线数据与预览验证

分块 NPY 与再次读取 RAW 的全部 840449 个 CD 事件逐事件一致。首 100ms
以 2ms 时间窗生成 50 帧 H.264 视频，1280×720、100fps，并通过视频重新解码
检查。输出为 `data/preview/G1_HARDWARE_CHECK_agent_raw_20261029_2110/`
中的 `preview_2ms.mp4` 和 `preview_2ms.png`。本次仅验证实际生成的 2ms 预览。

证据：`logs/agent_raw_validation_20261029_stats.json`、
`logs/agent_raw_validation_20261029_independent.json`、
`logs/agent_raw_validation_20261029_ffprobe.json`。
新增合法末尾纯时间字 EOF 回归后，全套 14 项测试通过，见
`logs/eof_regression.log`。

## 剩余验收

1. 相机专用静态地址已持久修好；在 Ubuntu 本地图形桌面点击 Device → OpenCamera，
   验证 Insight 实时显示、录制及原始 RAW 回放，当前等待用户反馈。
2. 排查 Python 实时 callback/batch-poll 超时；当前实时可用性仅由 C++ 路径确认。
3. 本轮等待 10 秒后无倒退不保证所有后续重启等待 10 秒均足够，旧传感器
   时钟状态可能不同。正式录制保留准备段并逐文件检查无倒退；
   填写镜头、光照、球体标记等实验元数据，
   完成真实旋转与 OptiTrack 真值比较。

初装时的合成验证事实保留在 `docs/verification_report.md`，不以本次硬件结果覆盖。
