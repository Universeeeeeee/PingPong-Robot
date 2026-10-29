"""Portable HTML report, generated only from measured benchmark JSON."""
from html import escape
from pathlib import Path


def write_demo(output, summaries):
    output = Path(output)
    render_animation(output)
    cards = []
    for s in summaries:
        name = escape(s['case'])
        rpm = '无效' if s['median_rpm'] is None else f'{s["median_rpm"]:.1f}'
        error = '—' if s.get('rpm_relative_error') is None else f'{100*s["rpm_relative_error"]:.2f}%'
        axis = '—' if s.get('median_axis_error_deg') is None else f'{s["median_axis_error_deg"]:.2f}°'
        cards.append(f'''<section id="{name}"><div class="heading"><h2>{name}</h2><a href="{name}/report.md">报告</a> · <a href="{name}/result.json">JSON</a> · <a href="{name}/windows.csv">CSV</a></div>
<div class="metrics"><div><b>{rpm}</b><span>RPM 中位数</span></div><div><b>{s['valid_fraction']:.0%}</b><span>有效窗口</span></div><div><b>{error}</b><span>RPM 相对误差</span></div><div><b>{axis}</b><span>有向轴误差</span></div></div>
<div class="plots"><img src="{name}/spherical_before_after.png" alt="球面补偿前后"><img src="{name}/rpm_curve.png" alt="转速曲线"></div><details><summary>优化过程与事件预览</summary><img src="{name}/optimization.png" alt="优化过程"><img src="{name}/event_preview.png" alt="事件预览"><p><a href="{name}/event_preview.mp4">打开 MP4</a></p></details></section>''')
    html = '''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>G1 球面旋转估计 · Demo</title>
<style>body{margin:0;background:#101a24;color:#e9f0f7;font:16px/1.6 system-ui,sans-serif}main{max-width:1120px;margin:auto;padding:40px 24px}h1{font-size:36px;margin:8px 0}h2{margin:0}a{color:#7ee2db}p{color:#b5c5d6}.badge{color:#7ee2db;letter-spacing:2px}.note{border-left:3px solid #f4b76b;padding:12px 20px;background:#25303a}section{margin:28px 0;padding:24px;background:#1b2a38;border-radius:16px}.heading{display:flex;gap:12px;align-items:center}.heading h2{flex:1}.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:16px;margin:24px 0}.metrics b{display:block;font-size:28px;color:#7ee2db}.metrics span{color:#b5c5d6}.plots{display:grid;grid-template-columns:1fr 1fr;gap:16px}img{max-width:100%;border-radius:8px;background:white}summary{cursor:pointer;margin:18px 0}details img{max-height:320px}nav{display:flex;flex-wrap:wrap;gap:16px;margin:24px 0}@media(max-width:700px){.plots,.metrics{grid-template-columns:1fr 1fr}.heading{display:block}h1{font-size:26px}}</style><main><div class="badge">DVSLume-G1 / PHASE 1 / OFFLINE</div><h1>固定球三维旋转估计</h1><p>球面提升 → SO(3) 反向补偿 → 球面 IWE → 多初值对比度优化</p><div class="note">当前展示理想球面特征轨迹的合成实验。真值仅用于评估；没有传给估计器。结果不能作为真实 G1 精度声明。真实 RAW 读取另行验证。</div>'''
    html += '<nav>' + ''.join(f'<a href="#{escape(s["case"])}">{escape(s["case"])}</a>' for s in summaries) + '</nav>'
    html += '<p>合成轨迹预览：1000 RPM，播放速度 0.1×（独立 300 ms 可视化片段）</p><img src="rotation.gif" alt="合成球旋转轨迹动画">'
    html += '<img src="benchmark.png" alt="各场景与真值对比">' + ''.join(cards)
    html += '<p>方法边界：工程多轴网格初始化；局部法向光流为单独对照。质量门限未经真实实验标定；低纹理、频闪、轴不可观测仍可能造成错误。</p></main></html>'
    (output/'index.html').write_text(html)


def render_animation(output):
    import cv2
    import numpy as np
    from PIL import Image
    from .synthetic import generate
    from .io import windows
    b = generate(rpm=1000, duration_s=.3)
    frames = []
    for start, end, batch in windows(iter([b]), 2000, start_us=0, duration_us=300000):
        frame = np.full((300, 300, 3), [24, 33, 45], dtype=np.uint8)
        cv2.circle(frame, (150, 150), 100, (85, 100, 115), 1)
        for polarity, color in [(0, (90, 180, 255)), (1, (255, 190, 70))]:
            for x, y in zip(batch.x[batch.polarity==polarity], batch.y[batch.polarity==polarity]):
                cv2.circle(frame, (round(x)+22, round(y)+22), 1, color, -1)
        cv2.putText(frame, f'SYNTHETIC 1000 RPM | {start/1000:.0f} ms', (9, 25), 0, .46, (220,230,240), 1)
        frames.append(Image.fromarray(cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)))
    frames[0].save(Path(output)/'rotation.gif', save_all=True, append_images=frames[1:], duration=20, loop=0)
