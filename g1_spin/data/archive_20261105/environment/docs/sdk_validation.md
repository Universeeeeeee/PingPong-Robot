# SDK 1.3.0 离线验证与边界

官方依据：

- https://sdk.dvsense.com/zh/html/install_zh.html
- https://sdk.dvsense.com/zh/html/file_zh.html
- https://sdk.dvsense.com/zh/html/python_dvs_sample_tutorial_zh.html
- https://github.com/dvsense/dvsense_driver_python_sample
- EVT3 测试编码：https://docs.prophesee.ai/stable/data/encoding_formats/evt3.html

主程序通过官方 Python `RawFileReader.load_file()`、`get_n_events()` 读取 RAW，
以官方 C++ `getMaxEvents()` 返回的索引事件总数作为停止条件，不自行解码相机文件。复制 SDK 返回的数组，
转换为无损整数结构 `t:uint64, x:uint16, y:uint16, p:uint8`。时间单位来自安装后的
`/usr/include/DvsenseBase/EventBase/EventTypes.hpp`：microseconds；极性 0 为亮度下降，1 为上升。
不排序、不减去起始时间、不修改坐标或极性。触发输入不是 CD 事件，本阶段不导出。

SDK 可能在 RAW 旁生成 `.dsmark.tmp` 等索引文件；不修改 RAW 本体。
SDK `load_file()` 自身进行全文件索引，耗时与文件大小有关；本工程的后续解码按批次保存。

## 已验证

`tests/make_evt3_fixture.py` 产生明确标注为 SYNTHETIC 的合法 EVT3 测试文件：
64×64、10000 个已知事件，时间戳从 123456 到 1123356 us，包含非零时间同步前导。
`test_official_raw_decoder_exact_values_and_chunk_tail` 对 250000/777 两种批大小执行
逐事件精确比较，并用独立编译、链接 apt SDK 的 C++ 程序复核数量及前几个事件。
时间同步前导不包含 CD 事件，所有 10000 个已知 CD 事件必须完整返回。

Python wheel 与 apt 库均为官方 1.3.0；Python 绑定在此版本中的回调名是
`add_event_stream_nocopy_callback`，不是由 C++ 方法名直接转成的名称。
`probe_camera.py` 已按安装后的绑定校正。2026-10-29 实机测试中，Python 实时
callback 与 batch-poll 路径均曾超时；C++ 实时取流通过，不能据此声称 Python 实时通过。

## 实测限制：人工 RAW 的索引起始块

在无时间同步前导的人工文件上，**load_file() 成功不代表解码完整**：

- 10000 个 CD 事件、约 80 KB 的文件报告 start=end=1123356 us，返回 0 个 CD 事件。
- 100000 个 CD 事件、每事件 10 us、首事件 123456 us 的文件，从 287296 us
  才开始返回事件。虽然返回数组长度仍为 100000，首 16384 个事件没有正确返回，
  数组还出现了 40 处时间戳倒退；仅核对数量会误判通过。本工程的时间戳检查会发现该问题。
- 同一现象可由 Python wheel 和链接系统 SDK 的 C++ 程序复现，不能归因于本工程的时间窗。
- 在人工流前加入合法的非零时间同步前导后，10000 个已知事件全部正确返回。

这些是对人工文件的实测结果，**不能据此断言真实 G1 固件录制文件一定丢失事件**，
也不能据此声称真实 RAW 已通过全量验收。工程保留相机 RAW 原件，且不会自动改写、
填充或伪造真实 RAW。索引总数尚未读满时若空读直接报错，避免静默产生不完整数据集。
首个真实录制应包含约 1 秒准备动作和随后实际实验，先通过 Insight 回放、事件统计、
SDK 索引起止时间和首次解码事件检查，再决定正式实验片段。准备动作不能代替完整性校验。
若真实文件出现同类问题，应保留原文件/校验值及最小复现交给厂商确认 SDK 行为。

## 二维 Baseline 的范围

Hough 圆候选（或人工测量圆）限定球体内部区域；Shi-Tomasi 特征 + 双向 LK 跟踪；
RANSAC 相似变换分离图像平移、均匀尺度与平面角度。角度除以实际事件时间窗中心间隔，
输出 rad/s 及有符号 rpm。图像 x 向右、y 向下，正号为从相机视角观察的顺时针。
每对至少 3 个内点；至少 5 个有效对才输出片段汇总值。失败窗口保留状态，不补造角速度。

合成数据是有纹理平面圆盘的事件采样，不是物理 DVS 球体仿真。已知真值单测验证
方向、单位、跟踪、噪声外点与视频导出流程，不能证明球面投影下的物理自旋精度。
完整三维自旋轴、镜头/几何标定和 OptiTrack 时间/坐标真值对齐留待下一阶段。

## 2026-10-29：真实 G1 RAW 的 EOF 与启动时钟

实际预热后录制文件 `G1_CONNECTIVITY_TEST_20261029_210504.raw` 的最后 CD 事件
为 12997713 us，SDK 文件结束时间为 12997718 us（尾部时间字）。在此文件上，
仅依赖 Python `reached_end_of_events()` 会一直为 false。Python 绑定未暴露
所需的 `getMaxEvents()`；`scripts/sdk_check --info RAW` 调用官方 C++ SDK 取得
索引总数，`event_io.py` 按剩余数量读取并在读满后停止。它不裁剪或重写 RAW，
索引为零、提前空读或返回数量超界均报错。

`logs/hardware_warmed_raw_stats.json` 记录 840449 个 CD 事件，与 C++ SDK 索引
总数一致；时间范围 9959239–12997713 us，时间戳倒退 0、非法事件 0。
这验证了该文件的离线读取与统计，不证明网络无丢包或物理旋转测量精度。
C++ 直播计数 860425 与录制文件统计覆盖区间不同，不作为逐事件等价比较。

首次立即录制的直播统计有 743 次倒退；独立 RAW 时间字检查发现 TIME_HIGH
从 1824 回到 0（约 7.47 秒回零），确认存在启动时钟重置。重新启动取流后
等待 10 秒再录 3 秒，本次直播与离线统计均未见倒退；这只是一次实测，
正式实验仍须逐文件检查。完整记录见 [硬件检查](connected_check_20261029.md)。

本次 EOF 回归使用含合法末尾纯时间字的合成 RAW：原 SDK 读完 10000 个
事件后 EOF 仍为 false；新读取路径以 250000/777 两种批大小均完整结束。
全套 14 项测试通过，日志为 `logs/eof_regression.log`。初装报告中的 13 项
为历史事实，保持不变。真实 RAW 分块 NPY 与再次读取 RAW 的 840449 个 CD
事件逐事件一致；首 100ms 的 2ms 预览生成 50 帧 H.264、1280×720、100fps，
并通过视频重新解码检查，见硬件检查报告中的证据路径。Insight GUI 回放仍待验证。
本轮等待 10 秒后无倒退，仅说明本轮测试成功；旧传感器时钟状态可能不同，
不能保证后续每次重启等待 10 秒均足够，正式录制须保留准备段并检查无倒退。
