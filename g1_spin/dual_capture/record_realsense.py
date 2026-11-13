"""Pre-roll RealSense recorder linked to Insight RAW creation, not hardware sync."""
import argparse,datetime,json,os,signal,struct,time,zlib
from pathlib import Path
import numpy as np
import pyrealsense2 as rs


def save(path,value):path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def png(path,array):
    a=np.asarray(array,dtype=np.uint8);h,w=a.shape
    def chunk(kind,data):return struct.pack('!I',len(data))+kind+data+struct.pack('!I',zlib.crc32(kind+data)&0xffffffff)
    path.write_bytes(b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('!2I5B',w,h,8,0,0,0,0))+chunk(b'IDAT',zlib.compress(b''.join(b'\0'+row.tobytes() for row in a)))+chunk(b'IEND',b''))

def raw_open(path):
    # Our user's Insight descriptors; returns unknown if its process cannot be inspected.
    found=False
    for proc in Path('/proc').glob('[0-9]*'):
        try:
            if 'DvsenseInsight' not in (proc/'comm').read_text():continue
            found=True
            for fd in (proc/'fd').iterdir():
                try:
                    if Path(os.readlink(fd))==path:return True
                except OSError:pass
        except (OSError,PermissionError):pass
    return False if found else None


def record(args):
    ctx=rs.context();devices=list(ctx.query_devices())
    if len(devices)!=1:raise RuntimeError(f'Expected exactly one RealSense, found {len(devices)}')
    device=devices[0];serial=device.get_info(rs.camera_info.serial_number)
    out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    watch=Path(args.watch).expanduser().resolve() if args.watch else None
    if watch and not watch.is_dir():raise ValueError('Insight output directory missing')
    known={str(p.resolve()) for p in watch.glob('*.raw')} if watch else set()
    meta=dict(status='starting',serial=serial,device=device.get_info(rs.camera_info.name),firmware=device.get_info(rs.camera_info.firmware_version),
              requested=dict(width=args.width,height=args.height,fps=args.fps,exposure_us=args.exposure_us),
              started_wall=datetime.datetime.now().astimezone().isoformat(),timestamp_unit_sdk='milliseconds',
              host_clock='monotonic_ns',sync='approximate_start_link_only; optical marker alignment required',
              bag='realsense.db3',watch_directory=str(watch),linked_raw=None)
    save(out/'session.json',meta);original=[];pipe=rs.pipeline();started=False;stopping=False;error=None;counts={};gaps={};last={};first_ts={};last_ts={};stats={};domains={};domain_changes={};first_host={};last_host={}
    def stop(*_):
        nonlocal stopping
        stopping=True
    signal.signal(signal.SIGINT,stop);signal.signal(signal.SIGTERM,stop)
    sensor=next(s for s in device.query_sensors() if s.supports(rs.option.emitter_enabled))
    options=[(rs.option.enable_auto_exposure,0),(rs.option.exposure,args.exposure_us),(rs.option.emitter_enabled,args.emitter),(rs.option.gain,args.gain)]
    try:
        meta['option_ranges']={}
        for opt,value in options:
            if not sensor.supports(opt):raise RuntimeError('Required sensor option missing: '+str(opt))
            rng=sensor.get_option_range(opt)
            meta['option_ranges'][str(opt)]={'min':rng.min,'max':rng.max,'step':rng.step}
            if not rng.min<=value<=rng.max:raise ValueError('Requested option out of range '+str(opt))
            original.append((opt,sensor.get_option(opt)));sensor.set_option(opt,value)
        meta['sensor_options']={str(opt):sensor.get_option(opt) for opt,_ in options}
        cfg=rs.config();cfg.enable_device(serial)
        cfg.enable_stream(rs.stream.infrared,1,args.width,args.height,rs.format.y8,args.fps)
        cfg.enable_stream(rs.stream.infrared,2,args.width,args.height,rs.format.y8,args.fps)
        cfg.enable_stream(rs.stream.depth,args.width,args.height,rs.format.z16,args.fps)
        cfg.enable_record_to_file(str(out/'realsense.db3'))
        meta['start_command_host_ns']=time.monotonic_ns();profile=pipe.start(cfg);started=True
        meta['start_return_host_ns']=time.monotonic_ns()
        profiles=[]
        for p in profile.get_streams():
            v=p.as_video_stream_profile();i=v.get_intrinsics();profiles.append(dict(stream=str(p.stream_type()),index=p.stream_index(),width=v.width(),height=v.height(),fps=p.fps(),format=str(p.format()),intrinsics=dict(fx=i.fx,fy=i.fy,ppx=i.ppx,ppy=i.ppy,model=str(i.model),coeffs=i.coeffs)))
        meta['profiles']=profiles
        if any(p['width']!=args.width or p['height']!=args.height or p['fps']!=args.fps for p in profiles):raise RuntimeError('Negotiated profile differs from request')
        streams=list(profile.get_streams());meta['extrinsics']=[]
        for a in streams:
            for b in streams:
                e=a.get_extrinsics_to(b);meta['extrinsics'].append(dict(source=[str(a.stream_type()),a.stream_index()],target=[str(b.stream_type()),b.stream_index()],rotation=e.rotation,translation=e.translation))
        meta['status']='recording';save(out/'session.json',meta)
        begin=time.monotonic();linked=None;linked_at=None;last_size=None;stable_since=None;log_at=begin;missing=0
        with (out/'frames.jsonl').open('w') as output,(out/'insight_events.jsonl').open('w') as links:
            while not stopping and time.monotonic()-begin<args.max_seconds:
                try:frames=pipe.wait_for_frames(2000);host=time.monotonic_ns();missing=0
                except RuntimeError:
                    missing+=1
                    if missing>=3:raise
                    continue
                for frame in frames:
                    p=frame.get_profile();key=f'{p.stream_type()}_{p.stream_index()}';n=frame.get_frame_number();ts=frame.get_timestamp()
                    if last.get(key)==n:continue
                    if key in last:gaps[key]=gaps.get(key,0)+max(0,n-last[key]-1)
                    current_domain=str(frame.get_frame_timestamp_domain())
                    if key in domains and current_domain!=domains[key]:domain_changes[key]=domain_changes.get(key,0)+1
                    first_host.setdefault(key,host);last_host[key]=host
                    last[key]=n;counts[key]=counts.get(key,0)+1;first_ts.setdefault(key,ts);last_ts[key]=ts;domains[key]=str(frame.get_frame_timestamp_domain())
                    md={}
                    for k in ['actual_exposure','gain_level','frame_timestamp','sensor_timestamp','backend_timestamp','time_of_arrival','actual_fps','auto_exposure']:
                        opt=getattr(rs.frame_metadata_value,k,None)
                        if opt is not None and frame.supports_frame_metadata(opt):md[k]=frame.get_frame_metadata(opt)
                    output.write(json.dumps(dict(stream=key,frame_number=n,sensor_timestamp_ms=ts,timestamp_domain=domains[key],host_monotonic_ns=host,metadata=md))+'\n')
                    if p.stream_type()==rs.stream.infrared and p.stream_index()==1 and counts[key] in [90,270,900,2700]:
                        image=np.asanyarray(frame.get_data());png(out/f'ir_left_{counts[key]:06d}.png',image)
                        stats[str(counts[key])]=dict(mean=float(image.mean()),p10=float(np.percentile(image,10)),p90=float(np.percentile(image,90)),saturated_fraction=float(np.mean(image>=250)))
                output.flush()
                now=time.monotonic()
                if watch:
                    new=sorted([p.resolve() for p in watch.glob('*.raw') if str(p.resolve()) not in known],key=lambda p:p.stat().st_mtime_ns)
                    if new and linked is None:
                        linked=new[0];linked_at=now;meta['linked_raw']=str(linked);meta['raw_seen_host_ns']=time.monotonic_ns();meta['raw_candidates']=[str(p) for p in new]
                        links.write(json.dumps(dict(kind='raw_created',path=str(linked),host_monotonic_ns=time.monotonic_ns()))+'\n');links.flush();save(out/'session.json',meta)
                        print('INSIGHT_RAW_DETECTED',linked,flush=True)
                    if linked:
                        size=linked.stat().st_size
                        if size!=last_size:stable_since=now;last_size=size
                        if stable_since is not None and now-stable_since>=3 and raw_open(linked) is False:
                            meta['stop_reason']='linked_raw_closed_and_stable';break
                        if now-linked_at>=args.max_raw_seconds:meta['stop_reason']='raw_duration_limit';break
                if now-log_at>=2:
                    print('RECORDING',round(now-begin,1),'s',counts,'raw=',str(linked),flush=True);log_at=now
            if stopping:meta['stop_reason']='signal'
            elif 'stop_reason' not in meta:meta['stop_reason']='duration_limit'
    except Exception as exc:
        error=repr(exc);meta['status']='failed';meta['error']=error
    finally:
        if started:
            pipe.stop();meta['stop_return_host_ns']=time.monotonic_ns()
        # Restore exposure before AE, which was saved first.
        restore_errors=[]
        for opt,value in reversed(original):
            try:sensor.set_option(opt,value)
            except Exception as exc:restore_errors.append(str(exc))
        meta.update(counts=counts,callback_frame_number_gaps=gaps,first_sensor_ms=first_ts,last_sensor_ms=last_ts,timestamp_domains=domains,
                    timestamp_domain_changes=domain_changes,observed_host_fps={k:(counts[k]-1)*1e9/(last_host[k]-first_host[k]) if last_host[k]>first_host[k] else None for k in counts},observed_fps={k:(counts[k]-1)*1000/(last_ts[k]-first_ts[k]) if last_ts[k]>first_ts[k] and not domain_changes.get(k) else None for k in counts},image_stats=stats,restore_errors=restore_errors)
        if not error:meta['status']='recorded_pending_alignment'
        save(out/'session.json',meta);print('DONE',out,meta['status'],flush=True)
    if error:raise RuntimeError(error)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);p.add_argument('--watch');p.add_argument('--max-seconds',type=float,default=120);p.add_argument('--max-raw-seconds',type=float,default=60)
    p.add_argument('--width',type=int,default=640);p.add_argument('--height',type=int,default=360);p.add_argument('--fps',type=int,default=90);p.add_argument('--exposure-us',type=float,default=500);p.add_argument('--emitter',type=int,choices=[0,1],default=0);p.add_argument('--gain',type=float,default=16);a=p.parse_args()
    if not 0<a.max_seconds<=300 or not 0<a.max_raw_seconds<=120:p.error('Recording duration outside bounded limits')
    record(a)
