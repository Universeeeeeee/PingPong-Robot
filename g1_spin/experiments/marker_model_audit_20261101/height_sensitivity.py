"""Analytic sensitivity of a raised material point; not a marker/DVS simulator."""
import csv
import hashlib
import json
import math
from pathlib import Path

OUT = Path(__file__).resolve().parent
BODY_RADIUS_MM = 20.0  # Nominal 40 mm diameter; not a calibrated measurement.
IMAGE_RADIUS_PX = 48.0  # Prior image estimate; not a verified body silhouette.
HEIGHTS_MM = [0.0, 0.1, 0.25, 0.5, 1.0]
IMAGE_RADII = [0.5, 0.8, 0.95]


def main():
    rows = []
    forward_error = 0.0
    for height in HEIGHTS_MM:
        scale = 1 + height / BODY_RADIUS_MM
        for rho in IMAGE_RADII:
            true_angle = math.asin(rho / scale)
            body_lift_angle = math.asin(rho)
            # Independent forward projection back to the chosen observed radius.
            forward_error = max(forward_error, abs(scale * math.sin(true_angle) - rho))
            rows.append(dict(
                radial_feature_height_mm=height,
                observed_radius_over_body_radius=rho,
                limb_projection_offset_bound_px=IMAGE_RADIUS_PX * height / BODY_RADIUS_MM,
                body_lift_direction_bias_deg=math.degrees(body_lift_angle - true_angle),
            ))
    assert forward_error < 1e-12
    assert all(abs(r['body_lift_direction_bias_deg']) < 1e-12 for r in rows if r['radial_feature_height_mm'] == 0)
    with (OUT / 'height_sensitivity.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = dict(
        model='fixed camera, orthographic projection, known body center, point at radius R+h',
        body_radius_mm=BODY_RADIUS_MM, image_radius_px=IMAGE_RADIUS_PX,
        heights_mm_are_hypothetical=True, forward_projection_max_error=forward_error,
        h_definition='radial height of observed material point above body, not marker diameter',
        exclusions=['real event generation', 'reflective photometry', 'marker silhouette',
                    'sphere occlusion', 'pose or spin estimation', 'measured height inference'],
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    (OUT / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print('| 假设径向高度 mm | 最大横向投影差 px | rho=0.5 提升方向偏差 deg | rho=0.8 | rho=0.95 |')
    print('|---:|---:|---:|---:|---:|')
    for height in HEIGHTS_MM:
        part = [r for r in rows if r['radial_feature_height_mm'] == height]
        print(f"| {height:g} | {part[0]['limb_projection_offset_bound_px']:.2f} | " +
              ' | '.join(f"{r['body_lift_direction_bias_deg']:.3f}" for r in part) + ' |')


if __name__ == '__main__':
    main()
