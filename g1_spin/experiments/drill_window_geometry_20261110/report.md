# 电钻球心与窗口离线对照结果

完成14个固定结束时刻 × 3种窗口 × 3种球心方案，共126条状态。生产算法未修改。

## 本轮结论

20ms + median3保留为下一轮候选。相对同档位名义参考，8个旋转片段中7个可用，平均绝对参考偏差11.52 RPM；10ms + fixed基线同7个可用片段为18.47 RPM。约38%的下降仅针对这组近似参考与已观察片段，不是独立真值精度提升认证。
evaluation块候选与基线各4/4可用，平均绝对参考偏差12.21与21.92 RPM；development块各3/4，为10.60与13.85 RPM。此处是结果分层，不声称提前锁定候选后盲测。
固定中心的20/40ms在3个有几何的纯平移窗口中均产生非零输出（37–67 RPM）；linear/median3在相同窗口没有接受非零输出。另1个平移窗几何失败，不能计作成功拒绝。静止两窗全部事件不足，不计作测准零。
40ms + median3仅6/8旋转窗口可用，参考偏差12.01 RPM，未体现整体优势；不直接把窗口越长当作更好。本轮仅优化既有几何观测的时间处理，未修复缺测的轮廓前端。

## 近似参考

用户已接受60fps手机标注作为近似转速参考。独立从marker CSV重算确认视频1为249.4359 RPM、视频2为244.3868 RPM。
用户本轮确认：另拍的同档位参考，转速近似相同。两段RPM平均246.9114作为名义档位参考，读数范围并集约233.88–261.53 RPM。该范围不是电钻跨录制漂移的标定界限。只报告相对档位参考偏差，不称同步真值误差；未把参考用于搜索、几何选择或调参。

## 分组结果

|录制/条件|窗口ms|球心|有效/预定|RPM中位数|有效范围|平均绝对参考偏差RPM|
|---|---:|---|---:|---:|---|---:|
|run1 rotation|10|fixed|4/4|249.4|222.5–265.3|18.3|
|run1 rotation|10|linear|4/4|273.2|228.6–285.9|27.5|
|run1 rotation|10|median3|4/4|261.7|228.6–266.0|16.7|
|run1 rotation|20|fixed|4/4|251.1|233.2–267.7|15.9|
|run1 rotation|20|linear|4/4|264.6|246.2–277.0|16.6|
|run1 rotation|20|median3|4/4|261.2|246.2–267.2|12.4|
|run1 rotation|40|fixed|3/4|255.2|237.2–262.8|11.3|
|run1 rotation|40|linear|3/4|262.8|255.7–271.3|16.4|
|run1 rotation|40|median3|3/4|261.6|255.7–269.8|15.5|
|run4 rotation|10|fixed|3/4|258.3|235.6–280.5|18.7|
|run4 rotation|10|linear|3/4|253.7|221.7–254.8|13.3|
|run4 rotation|10|median3|3/4|253.7|229.4–255.2|10.9|
|run4 rotation|20|fixed|3/4|259.5|228.7–275.1|19.7|
|run4 rotation|20|linear|3/4|244.2|223.1–257.7|12.4|
|run4 rotation|20|median3|3/4|243.6|229.7–257.6|10.4|
|run4 rotation|40|fixed|3/4|271.2|229.0–275.6|23.6|
|run4 rotation|40|linear|3/4|252.7|248.9–265.7|8.9|
|run4 rotation|40|median3|3/4|252.3|243.2–263.4|8.6|
|run5 translation|10|fixed|0/4|—|—|—|
|run5 translation|10|linear|0/4|—|—|—|
|run5 translation|10|median3|0/4|—|—|—|
|run5 translation|20|fixed|3/4|60.8|45.2–67.2|—|
|run5 translation|20|linear|0/4|—|—|—|
|run5 translation|20|median3|0/4|—|—|—|
|run5 translation|40|fixed|3/4|59.4|37.4–61.0|—|
|run5 translation|40|linear|0/4|—|—|—|
|run5 translation|40|median3|0/4|—|—|—|
|run5 static|10|fixed|0/2|—|—|—|
|run5 static|10|linear|0/2|—|—|—|
|run5 static|10|median3|0/2|—|—|—|
|run5 static|20|fixed|0/2|—|—|—|
|run5 static|20|linear|0/2|—|—|—|
|run5 static|20|median3|0/2|—|—|—|
|run5 static|40|fixed|0/2|—|—|—|
|run5 static|40|linear|0/2|—|—|—|
|run5 static|40|median3|0/2|—|—|—|

## 计算与解释边界

单次有输入估计耗时中位数2.289s，最大7.844s；包含球面提升和搜索，不含文件读取/原轮廓提取。
fixed：窗口中点固定中心；linear：逐事件插值；median3：连续三点中值平滑后插值。三者同窗共享事件ID，40ms等长窗会因几何缺口失败，失败保留。
各长度采用相同窗口结束时刻，但事件数和5000上限造成采样密度差异；结果是当前计算预算下整体表现，不能单独归因于时长。
几何采用离线观测，可能含未来支点；lookahead字段是三方案共同的保守上界，非在线延迟。半径保持60px，未优化半径、外轮廓检测器或CMAX门限。
development/evaluation按既有时间块区分，本次方案预先固定；历史已观察过这些片段，不能冒称独立盲测。
未输出不等于测到零；转速稳定不等于轴正确。三维轴仅保留候选，未报告精确轴误差。

## 复现

`g1_spin/.venv/bin/python g1_spin/experiments/drill_window_geometry_20261110/run.py`
首次运行输出拒绝覆盖。后续执行同目录 verify.py 与 summarize.py 核验并生成汇总。
run.py保存运行前未配对的参考状态；reference.json另记录运行期间用户确认的同档位关系及参考汇总。重新开展实验时需显式保留该确认，summarize.py要求确认后的名义参考字段。
来源与参数见protocol.md/freeze.json，逐条结果见results.csv/results.jsonl，独立输入核验见independent_verification.json。
