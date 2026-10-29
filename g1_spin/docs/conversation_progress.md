# 本轮进展归档（2026-10-29）

## 实施与交付

1. 检查 Ubuntu 22.04 / Python 3.10 / 官方 DVSense SDK 1.3.0，确认已有工程和真实连通性RAW。
2. 在 Ubuntu `~/g1_spin/phase1_demo` 独立开发，原有SDK读取程序及数据未覆盖。
3. 实现事件格式/读取桥接、球面几何、周期Baseline、局部法向光流、s-CMax、多初值质量检查、CLI、合成器和报告。
4. 三个子agent并行处理算法、数据/测试、文献核对；核对Sony公开伪代码笔误并建立资料映射。
5. 合成demo及RAW读取验收见 `acceptance.md`；原工程14项回归通过。
6. Ubuntu从172.16.16.225改用192.168.119.157连接后，将三段新录制RAW和官方SDK导出的908个NPY块复制到Mac，逐文件SHA256验证。
7. Mac独立Python3.12环境支持manifest分块读取，无需远程连接。最新测试30项通过。
8. 最新17.58秒录像选取1.6–1.9、6.2–6.5、11.2–11.5秒进行分析，得到约192、137、76 RPM候选；独立图案周期结果接近。偏心运动、遮挡和ROI不确定性使三维轴尚未获得物理验证。

## 仓库内交付入口

- `../README.md`：安装与运行。
- `../outputs/demo_final/index.html`：合成demo。
- `../outputs/real_222403/index.html`：真实录像分析，包含实际事件GIF、候选减速曲线、补偿图与诊断。
- `../outputs/older_recordings_qc/report.md`：另外两段录像筛选。
- `../data/archive_manifest.json`：三段RAW路径、长度、事件数和SHA256。
- `../data/decoded/*/manifest.json`、SHA256SUMS与校验记录：SDK解码来源及完整性证据。
- `phase1_task.md`、`reference_brief.md`：用户原始任务与补充资料。

## 未提交的大数据与复现条件

仓库现有规则排除大规模数据及虚拟环境。约3.6GB RAW/NPY仍保存在本机：
`/Users/vae/Projects/PingPongRobot/g1_spin/data`，未删除。Ubuntu原件仍在`~/g1_spin/data/raw`。
本次仅提交数据元信息及哈希，不将大型RAW/NPY塞入普通Git。
恢复真实分析需要从上述本地备份或Ubuntu取得原始数据块，并对照保存的校验清单；仓库克隆本身不包含原始录像。
`outputs/demo_final` 中小规模合成NPZ、图表、MP4及报告作为已执行demo的证据纳入。
未保存到Mac的1.85GB录像（221635）不在三段归档清单中。

## 当前边界

完整来源、实验条件与风险分别见`references.md`、`experiment_protocol.md`和真实分析报告。
s-CMax采用工程多轴网格初值，光流是独立对照，尚未连接为初始化。合成数据为理想轨迹而非物理DVS模拟。
真实录像没有独立转速/轴真值；不能据优化收敛声称可信三维轴或论文精度，也不能把末尾事件减少当作0RPM。
本次提交仅覆盖g1_spin目录，不包含仓库中原有无关README修改、Optitrack/mocap等删除或其他工程文件。
