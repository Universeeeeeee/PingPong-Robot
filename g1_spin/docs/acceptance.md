# Demo 验收记录 · 2026-10-29

## 交付位置与环境

Ubuntu `pingpong1@172.16.16.225:~/g1_spin/phase1_demo`，独立 `.venv`。
Ubuntu 22.04.5、Intel Core Ultra 9 285HX（24 CPU）、30 GiB 内存、约 870 GiB 可用磁盘，Python 3.10.12。
原工程已安装官方 DVSense SDK/Insight 1.3.0；本次复用其 RAW 读取适配器，没有变更网络、固件、原程序或 RAW 内容。
最终环境不继承系统 site-packages；`pip check` 通过，版本固定在 `requirements-lock.txt`。
初次 PyPI 下载超时后使用清华镜像完成安装；安装阶段 setuptools/packaging 不兼容已在独立环境修复。

## 实际执行

```bash
cd ~/g1_spin/phase1_demo
.venv/bin/python -m pytest -q
# 29 passed in 4.76s
.venv/bin/python -m pip check
# No broken requirements found.
.venv/bin/python -m g1_spin.cli synthetic-benchmark --output outputs/demo_final
.venv/bin/python -m g1_spin.cli inspect ../data/raw/G1_CONNECTIVITY_TEST_20261029_210504.raw --output outputs/raw_check
.venv/bin/python -m g1_spin.cli preview ../data/raw/G1_CONNECTIVITY_TEST_20261029_210504.raw --output outputs/raw_check/preview.mp4 --duration-ms 100
cd ~/g1_spin
.venv/bin/python tests/test_offline.py
# 14 tests, OK（原工程回归）
```

新增测试覆盖：整数时间差、微秒单位、非法输入、球面提升、SO(3) 正逆方向、phi 接缝与 IWE 质量守恒、已知转速恢复、静止/噪声/退化拒绝、周期信号、滑动窗、空窗、光流精确约束与秩退化、实际 MP4 解码及报告。
这不是 P0–P8 全项完成声明。

## 合成实测（10 ms，3 窗/场景）

| 场景 | 有效窗口 | RPM 中位数 | RPM 相对误差 | 轴误差中位数 |
|---|---:|---:|---:|---:|
| 500 RPM | 3/3 | 501.015 | 0.203% | 4.182° |
| 1000 RPM | 3/3 | 1022.214 | 2.221% | 1.739° |
| 2000 RPM | 3/3 | 2010.344 | 0.517% | 0.237° |
| 3000 RPM | 3/3 | 3000.034 | 0.0011% | 0.850° |
| 反向 1000 RPM | 3/3 | 997.159 | 0.284% | 0.561° |
| 10% 背景噪声 | 3/3 | 1022.184 | 2.218% | 1.738° |
| 空输入 | 0/3 | null | — | — |

每场景优化总耗时约 0.8–0.9 s（不含绘图/视频），不是实时性能声明。500 RPM 的一个窗口轴误差为 5.485°，没有达到逐窗 <5° 调试目标；表中是中位数，不能据此掩盖该窗口。

初版 `outputs/demo` 使用 4200 Hz 合成采样，3000 RPM 可见轨迹较少，每窗仅84个有效事件，被100事件质量门限拒绝。保留该输出；最终所有场景统一改为8400 Hz合成采样，没有降低门限。采样密度是理想生成器工程参数，不是 G1 实际响应模型。

周期法在独立0.6秒样本上得到500/1000/2000/反向/噪声候选频率；3000 RPM轴向导致无可靠事件率调制而无效。周期候选不等于已证明的一圈频率。光流对照在全部点轨迹场景缺少稳定二维局部平面，均返回无效；数学约束测试通过不等于真实事件流精度验证。三种方法目前并非同等时长的完整性能对照，数据见 `benchmark.json`。

## RAW 与可视化

- 完整读取840449个CD事件；1280×720；9959239–12997713 us，3.038474 s。
- 时间戳倒退0；相同时间戳107022对属合法现象；最大间隔49 us。
- RAW前后SHA256相同，记录在 `outputs/raw_check/raw_sha256_{before,after}.txt`。
- 7份合成MP4均逐帧重新解码为15帧256×256；真实RAW预览重新解码为50帧1280×720。证据 `outputs/video_validation.json`。
- HTML中另有150帧、3秒循环GIF，展示独立300ms合成片段的0.1倍速播放。
- 1000 RPM补偿前后图已人工查看：轨迹集中为点；两幅采用相同色标。

现有RAW是连通性录制，没有旋转真值、球ROI或静止球实验标签，未对其虚构球转速。

## 剩余工作与风险

1. 正式固定球录制与独立转速/轴真值；真实旋转精度尚未验证。
2. 光流前端在实际商标边缘上验证，并作为s-CMax初值；当前核心使用工程多轴搜索。
3. 更真实的对数亮度事件模拟、局部商标不可见/遮挡/频闪/热点负例；当前全局离散标记轨迹容易于真实商标。
4. 局部可观测性和质量门限需要校准；当前多初值、空间分布、相邻窗检查不能保证轴可辨识。
5. GUI点选ROI尚未实现；当前直接输入球心/半径，提供双圆预览。热点检测和背景干扰自动拒绝未实现，outside_roi_fraction仅为统计。
6. 尚未形成同一有效事件序列、同等观测时长上的A/B/C完整算法比较。
7. 当前64×128经纬网格存在极区畸变；搜索上限可配置，最小转速、网格与质量门限尚未全部配置化。
8. RAW桥接依赖原工程官方SDK辅助程序；NPZ是小型演示格式，整份加载，不能冒称大文件流式格式。

## 下一步只需用户做的事

打开 `outputs/demo_final/index.html` 演示；随后按 `experiment_protocol.md` 录制有非对称商标的固定球，提供独立转速参考、球心半径及实验元数据。配置好ROI后即可用统一 `analyze` 命令处理，不需要改算法代码。
