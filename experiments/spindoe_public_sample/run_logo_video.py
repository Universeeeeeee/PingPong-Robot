"""Run print-feature tracking on a video with explicit ball circles.

Without --capture-fps, report radians per VIDEO FRAME, never physical rad/s.
The container playback FPS is deliberately not used as the measurement clock.
"""
import argparse
import hashlib
import json
from pathlib import Path

import cv2
import imageio.v2 as imageio
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial.transform import Rotation

from logo_frontend import CENTER, RADIUS, SIZE, _normalized_frame, estimate_sequence


def run(video, circles_path, output, start=0, count=120, capture_fps=None):
    if capture_fps is not None and (not np.isfinite(capture_fps) or capture_fps <= 0):
        raise ValueError("capture_fps must be positive and refer to acquisition")
    cap = cv2.VideoCapture(str(video))
    if not cap.isOpened():
        raise ValueError(f"Cannot decode {video}")
    container_fps = cap.get(cv2.CAP_PROP_FPS)
    cap.set(cv2.CAP_PROP_POS_FRAMES, start)
    frames = []
    for _ in range(count):
        ok, frame = cap.read()
        if not ok:
            break
        frames.append(cv2.cvtColor(frame, cv2.COLOR_BGR2RGB))
    cap.release()
    if len(frames) < 2:
        raise ValueError("Selected video interval has fewer than two decoded frames")
    geometry = json.loads(circles_path.read_text())
    if "circle" in geometry:
        circles = np.tile(geometry["circle"], (len(frames), 1))
    else:
        circles = np.asarray(geometry["circles"])
        if geometry.get("start_frame", start) != start or len(circles) != len(frames):
            raise ValueError("Circle annotations must match the selected frame interval")
    intervals = 1/capture_fps if capture_fps else 1.0
    timestamps = np.arange(len(frames))*intervals
    result = estimate_sequence(frames, timestamps, circles)
    result.update(source_video=str(video.resolve()), source_sha256=hashlib.sha256(video.read_bytes()).hexdigest(),
                  first_source_frame=start, decoded_frame_count=len(frames), container_playback_fps=container_fps,
                  supplied_capture_fps=capture_fps,
                  time_basis="supplied_capture_fps" if capture_fps else "video_frame_index",
                  geometry_annotation=geometry)
    result["preview_fps"] = 2 if geometry.get("kind") == "real image strip replay" else 20
    for window in result["windows"]:
        omega = window["omega_rad_s"]
        window["step_angle_deg"] = None if omega is None else float(np.degrees(np.linalg.norm(omega)*intervals))
        if capture_fps is None:
            window["rotation_vector_rad_per_video_frame"] = window.pop("omega_rad_s")
            window.pop("start_time_s")
            window.pop("end_time_s")
    output.mkdir(parents=True, exist_ok=True)
    (output/"results.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    font = ImageFont.load_default(size=17)
    small = ImageFont.load_default(size=14)
    rendered=[]
    for i, frame in enumerate(frames):
        normalized = _normalized_frame(frame, circles[i])
        raw = cv2.cvtColor(normalized, cv2.COLOR_GRAY2RGB)
        tracked, predicted = raw.copy(), raw.copy()
        window = result["windows"][i-1] if i else None
        if window:
            a = np.asarray(window["start_points_px"],dtype=float).reshape(-1,2)
            b = np.asarray(window["end_points_px"],dtype=float).reshape(-1,2)
            a = (a-circles[i-1,:2])*RADIUS/circles[i-1,2]+CENTER
            b = (b-circles[i,:2])*RADIUS/circles[i,2]+CENTER
            for previous, current in zip(a,b):
                cv2.arrowedLine(tracked,tuple(np.rint(previous).astype(int)),tuple(np.rint(current).astype(int)),(40,230,240),1,tipLength=0.25)
                cv2.circle(predicted,tuple(np.rint(current).astype(int)),1,(60,240,80),-1)
            vector = window.get("omega_rad_s",window.get("rotation_vector_rad_per_video_frame"))
            if vector is not None and len(a):
                xy=(a-CENTER)/RADIUS
                xy[:,1]*=-1
                points=np.c_[xy,np.sqrt(np.maximum(0,1-np.sum(xy*xy,axis=1)))]
                moved=Rotation.from_rotvec(np.asarray(vector)*intervals).apply(points)
                projected=np.c_[moved[:,0],-moved[:,1]]*RADIUS+CENTER
                for point,keep in zip(projected,window["inlier_mask"]):
                    if keep:
                        cv2.circle(predicted,tuple(np.rint(point).astype(int)),1,(255,80,160),-1)
        canvas=Image.new("RGB",(960,420),(18,22,29))
        draw=ImageDraw.Draw(canvas)
        clock="capture FPS supplied" if capture_fps else "capture timing UNKNOWN; no rps"
        kind=geometry.get("kind","video input").upper()
        draw.text((20,12),f"{kind} | print features -> sphere motion | {clock}",font=small,fill="white")
        for column,(label,panel) in enumerate([("Normalized input",raw),("Tracked print features",tracked),("Measured green / fitted pink",predicted)]):
            x=column*320+22
            draw.text((x,52),label,font=small,fill="white")
            canvas.paste(Image.fromarray(panel).resize((272,272),Image.Resampling.NEAREST),(x,80))
        text=f"Frame {start+i} | initialization" if window is None else f"Frame {start+i} | {window['reason']} | tracks {window['track_count']} | inliers {window['inlier_count']}"
        draw.text((22,370),text,font=font,fill="white")
        detail="No rotation estimate" if not window or not window["valid"] else f"Step: {window['step_angle_deg']:.3f} deg / video frame; fit RMS: {window['reprojection_rms_px']:.3f} normalized px"
        draw.text((22,398),detail,font=small,fill=(180,190,205))
        rendered.append(canvas)
    # Display cadence is independent of the measurement time base.
    with imageio.get_writer(str(output/"tracks_preview.mp4"),fps=result["preview_fps"],codec="libx264",quality=8,macro_block_size=4) as writer:
        for frame in rendered:
            writer.append_data(np.asarray(frame))
    selected=min(len(frames)//2,len(frames)-1)
    rendered[selected].save(output/"preview.png")
    Image.fromarray(frames[selected]).save(output/"source_frame.png")
    print(json.dumps({k:result[k] for k in ["status","valid_window_count","window_count","time_basis","decoded_frame_count"]},indent=2))
    return result


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("video",type=Path)
    parser.add_argument("--circles",type=Path,required=True,help="JSON: circle [cx,cy,r] or circles Nx3 plus start_frame")
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--start",type=int,default=0)
    parser.add_argument("--count",type=int,default=120)
    parser.add_argument("--capture-fps",type=float,help="Explicit acquisition cadence; never infer from MP4 playback")
    args=parser.parse_args()
    run(args.video,args.circles,args.out,args.start,args.count,args.capture_fps)
