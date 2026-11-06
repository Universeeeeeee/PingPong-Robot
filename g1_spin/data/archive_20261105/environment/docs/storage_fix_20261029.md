# Insight 保存目录修复（2026-10-29）

实际保存目录：`/home/pingpong1/g1_spin/data/raw`。视频导出：`/home/pingpong1/g1_spin/data/preview`。已备份并修改 config.ini，然后重启 Insight。启动后的保存路径错误消失，用户完成 GUI 录制后已验证真实文件：

- `DSP03F100025_20261029_221544.raw`
- 大小：716230574 字节。
- 官方 SDK 读取 186030505 个事件。
- 事件时间跨度约 22.654 秒，1280×720，时间戳倒退 0，非法事件 0。
- 桌面目录入口：`/home/pingpong1/Desktop/G1录制文件`。

旧配置曾指向不存在的 `/opt/DvsenseInsight/Data` 及 `/home/DvdenseInsight/Data`。在 /home、/opt/DvsenseInsight、/tmp、/media、/mnt 的可访问范围未找到此前22:03、22:04、22:06的片段；没有可访问的已删除但仍打开的 RAW。未声称这些片段已恢复。

证据：`logs/storage_fix.json`、`logs/insight_saved_recording_stats.json`、`logs/insight_storage_fixed.log`。GUI 录制落盘已验证；GUI 打开 RAW 回放不在本次验证范围。
