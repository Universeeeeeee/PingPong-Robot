# B0 已有 Marker 布局与标定资产核查

日期：2026-11-01。仅检查本地工作区；未连接Motive、动捕电脑或设备，未改资产。

## 发现了布局文件，但尚未绑定当前球

文件：`HOPE/mocap/Optitrack Ping Pong Ball Marker Template.csv`。

- 格式首行：`File Type, RigidBody Asset`，Version 1。
- 名称：`Ball_40mm_10pt`。
- 单位：`m`，StreamingID 24。
- 包含Point 1至Point 10的10组三维坐标，各点到资产原点的距离约20mm。
- 点集非共面（中心化坐标矩阵秩3）；文件描述三维点阵，不只是一个球心位置。
- HOPE子仓库记录：`6a39688`，2026-08-03创建该文件。
- 本次文件SHA256：`73cb42a43d65454561474946bac2ef05c86b6f707a4b6b30899a549d89e3d724`。

用户报告2026-10-29的G1录像使用约14个凸起球形反光Marker。没有找到把该10点文件与这颗实物球、当前Marker编号、实际几何中心或10月29日录制绑定的记录。因此它是可读的参考模板，不是本次实验可直接采用的已知三维布局。约14点是用户估计，也不据此断言两者绝不可能有关；需要实物/资产记录确认。

## 其他证据及其边界

- `HOPE/docs/OPTITRACK.md`和`HOPE/mocap/HOPE_Motion_Capture_System_and_Coordinates_Reference_Setup.md`描述了名为Ball的6-DOF刚体、球心pivot与非对称点阵。这是系统配置/规范，不是当前G1录制的标定结果。
- `7.19-7.25/windows部署/optitrack_windows_bundle_20260728/WINDOWS_QUICKSTART.md`记录历史现场曾使用ball/table刚体。它没有提供当前约14点的局部三维坐标。
- 本地工作区按扩展名检查（包括隐藏及忽略文件，排除.git、虚拟环境、node_modules、第三方目录），未找到`.motive/.tak/.tra/.mdef/.cal/.ttp`资产或录制文件。名称/内容检索也未找到与该G1录制绑定的约14点布局。
- 在`g1_spin/configs`、`g1_spin/docs`、`g1_spin/data`的JSON/YAML/Markdown中，未发现可绑定本次录制的相机内参矩阵和畸变参数。其他相机或机器人标定不自动适用于G1。

这些是限定本地范围内的检索结果，不代表Motive电脑或其他未搜索存储中没有资产。

## 后续能复用的条件

若取得本次实际Ball资产的导出，应核对：Marker数量与ID、刚体局部三维坐标、单位、pivot与轴向、是否为Marker几何中心、实物布置是否改变，以及G1相机内参/畸变。几何布局可用于建立后续重投影模型；被核验系统的旋转输出不能反过来作为事件估计的独立正确性证明。

B0的二维局部观测与短时关联不依赖这些资产，可以继续；当前不满足直接进入已知布局PnP的输入条件。
