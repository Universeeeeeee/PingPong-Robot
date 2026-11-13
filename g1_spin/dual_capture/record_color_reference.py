"""Bounded color reference recording, projector off; Insight owns event camera."""
import sys,time,json,signal,datetime
from pathlib import Path
import numpy as np
import pyrealsense2 as rs
from record_realsense import png,raw_open,save
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
watch=Path.home()/'g1_spin/data/raw';known={str(p.resolve()) for p in watch.glob('*.raw')}
p=rs.pipeline();cfg=rs.config();cfg.enable_stream(rs.stream.color,640,360,rs.format.rgb8,90);cfg.enable_record_to_file(str(out/'color.db3'))
stop=False

def halt(*_):
 global stop
 stop=True
signal.signal(signal.SIGINT,halt);signal.signal(signal.SIGTERM,halt)
profile=None;meta={'status':'starting','projector_requested':0,'stream':'color 640x360 90fps','clock_alignment':'unverified optical markers required'}
try:
 devs=list(rs.context().query_devices())
 if len(devs)!=1:raise RuntimeError('Expected one RealSense')
 devs[0].first_depth_sensor().set_option(rs.option.emitter_enabled,0)
 profile=p.start(cfg);sensor=profile.get_device().first_color_sensor()
 v=profile.get_stream(rs.stream.color).as_video_stream_profile()
 if (v.width(),v.height(),v.fps())!=(640,360,90):raise RuntimeError('Profile mismatch')
 meta['actual_profile']={'width':v.width(),'height':v.height(),'fps':v.fps(),'format':str(v.format())}
 sensor.set_option(rs.option.enable_auto_exposure,1)
 for _ in range(60):f=p.wait_for_frames(2000).get_color_frame()
 exposure=min(sensor.get_option(rs.option.exposure),50.0);gain=sensor.get_option(rs.option.gain)
 sensor.set_option(rs.option.enable_auto_exposure,0);sensor.set_option(rs.option.exposure,exposure);sensor.set_option(rs.option.gain,gain)
 sensor_times=[];frame_nums=[];domains=[];host0=time.monotonic()
 for _ in range(180):
  f=p.wait_for_frames(2000).get_color_frame();sensor_times.append(f.get_timestamp());frame_nums.append(f.get_frame_number());domains.append(str(f.get_frame_timestamp_domain()))
 elapsed=time.monotonic()-host0
 dt=np.diff(np.array(sensor_times)[-90:]);measured=1000/float(np.median(dt));host_fps=180/elapsed
 meta['preflight']={'median_frame_ms':float(np.median(dt)),'sensor_fps':measured,'host_fps':host_fps,'frame_gaps':int(np.maximum(np.diff(frame_nums)-1,0).sum()),'timestamp_domains':sorted(set(domains))}
 if not 85<=measured<=95 or host_fps<80:raise RuntimeError('90FPS delivery check failed: '+str(meta['preflight']))
 meta.update(status='recording',exposure_sdk_units=exposure,gain=gain,auto_exposure=0,wall_start=datetime.datetime.now().astimezone().isoformat(),serial=profile.get_device().get_info(rs.camera_info.serial_number))
 save(out/'session.json',meta)
 im=np.asanyarray(f.get_data());png(out/'preview_gray.png',im.mean(axis=2).astype(np.uint8));print('READY',out,'profile',meta['actual_profile'],'preflight',meta['preflight'],'exposure_sdk',exposure,'mean',float(im.mean()),flush=True)
 begin=time.monotonic();linked=None;lastsize=None;changed=begin;linkedtime=None;count=0;gaps=0;last=None;log=begin
 with (out/'frames.jsonl').open('w') as stream:
  while not stop and time.monotonic()-begin<300:
   f=p.wait_for_frames(3000).get_color_frame();now=time.monotonic();n=f.get_frame_number()
   if last==n:continue
   if last is not None:gaps+=max(0,n-last-1)
   last=n;count+=1;md={}
   for key in ['sensor_timestamp','frame_timestamp','actual_exposure','gain_level','auto_exposure']:
    v=getattr(rs.frame_metadata_value,key,None)
    if v is not None and f.supports_frame_metadata(v):md[key]=f.get_frame_metadata(v)
   stream.write(json.dumps(dict(frame_number=n,sdk_ms=f.get_timestamp(),domain=str(f.get_frame_timestamp_domain()),host_ns=time.monotonic_ns(),metadata=md))+'\n');stream.flush()
   new=[q.resolve() for q in watch.glob('*.raw') if str(q.resolve()) not in known]
   if linked is None and new:
    linked=sorted(new,key=lambda q:q.stat().st_mtime_ns)[0];linkedtime=now;meta.update(linked_raw=str(linked),seen_host_ns=time.monotonic_ns());save(out/'session.json',meta);print('LINKED',linked,flush=True)
   if linked:
    size=linked.stat().st_size
    if size!=lastsize:lastsize=size;changed=now
    if now-changed>3 and raw_open(linked) is False:meta['stop_reason']='raw_closed';break
    if now-linkedtime>120:meta['stop_reason']='120s_event_limit';break
   if now-log>5:print('RECORDING',count,'linked',str(linked),flush=True);log=now
 meta.update(status='recorded_pending_verification',color_received=count,callback_frame_gaps=gaps)
except Exception as e:
 meta.update(status='failed',error=repr(e));print('ERROR',repr(e),flush=True)
finally:
 if profile:p.stop()
 save(out/'session.json',meta);print('DONE',meta,flush=True)
