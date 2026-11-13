"""Phase annotation and coarse paired review; no precise synchronization claim."""
import json,sqlite3,hashlib
from pathlib import Path
import numpy as np,cv2
from audit_remote import OUT,SESSION,RAW,parse_image,save

# User reports order; numerical boundaries are visual, approximate and not angular truth.
PHASES=[(0,3,'shared hand marker / preparation','visual'),(3,7.5,'horizontal translation','user order + visual'),(7.5,13,'vertical forward rotation','user order + shaft pose'),(13,14.5,'pose transition','visual'),(14.5,21.5,'axis toward camera','user order + visual'),(21.5,23,'pose transition','visual'),(23,29,'tilted forward rotation','user order; reversal boundary unverified'),(29,35,'tilted reverse rotation','user order; reversal boundary unverified'),(35,36.5,'pose transition','visual'),(36.5,42,'vertical reverse rotation','user order + shaft pose'),(42,52,'multi-direction shaking','user order + visual'),(52,54.064668,'end movement','visual')]

def main():
    offset=json.loads((OUT/'rough_alignment.json').read_text())['corresponding_db_s'];origin=json.loads((OUT/'event_stats.json').read_text())['first_timestamp_us']
    phases=[dict(event_start_s=a,event_end_s=b,d455_db_start_s=a+offset,d455_db_end_s=b+offset,label=label,evidence=evidence,boundary_precision='coarse visual, roughly 0.5-1s; tilted reversal time unresolved',spin_ground_truth=False) for a,b,label,evidence in PHASES]
    save('phase_annotations.json',dict(user_description='水平移动几秒→原地竖直转→朝相机转→斜转→斜向反转→竖直反转→各方向摇晃→结束',time_mapping='d455_db_seconds ~= event_seconds_since_first_CD + 205.443818693',mapping_kind='file_seen_host_anchor_plus_visual_consistency; NOT calibrated sync',offset_s=offset,drift_estimated=False,phases=phases))
    metadata=json.loads((OUT/'ir_left_metadata.json').read_text());active=[r for r in metadata if offset<=r['db_ns']/1e9<=offset+54.064668]
    n=np.array([int(r['Frame number']) for r in active]);t=np.array([int(r['Sensor Timestamp']) for r in active]);dt=np.diff(t)
    audit=dict(frames=len(active),frame_numbers_first_last=[int(n[0]),int(n[-1])],missing_frame_numbers=int(np.maximum(np.diff(n)-1,0).sum()),sensor_dt_us_quantiles=np.percentile(dt,[0,50,90,99,100]).tolist(),timestamp_domains=list(dict.fromkeys(r['timestamp_domain'] for r in active)),actual_exposure_us=sorted(set(int(r['Actual Exposure']) for r in active)),effective_sensor_hz=(len(t)-1)*1e6/(t[-1]-t[0]),recording_span_s=(t[-1]-t[0])*1e-6)
    save('active_interval_d455_audit.json',audit)
    anomalies=json.loads((OUT/'event_timestamp_anomalies.json').read_text());intervals=sorted([((min(r['previous'],r['current'])-origin)/1e6-.02,(max(r['previous'],r['current'])-origin)/1e6+.02) for r in anomalies]);merged=[]
    for a,b in intervals:
        a=max(0,a);b=min(54.064668,b)
        if merged and a<=merged[-1][1]:merged[-1][1]=max(merged[-1][1],b)
        else:merged.append([a,b])
    save('event_exclusion_intervals.json',dict(reason='20ms padding around decoded timestamp regressions; diagnostic only, no timestamp repair',intervals_s=merged,total_excluded_s=sum(b-a for a,b in merged),timestamp_regressions=50,max_backward_us=max(r['previous']-r['current'] for r in anomalies)))
    db=sqlite3.connect('file:'+str(SESSION/'realsense.db3')+'?mode=ro',uri=True)
    writer=cv2.VideoWriter(str(OUT/'paired_coarse_review.mp4'),cv2.VideoWriter_fourcc(*'mp4v'),2,(1280,420));assert writer.isOpened();contact=[]
    for path in sorted((OUT/'event_samples').glob('*.png')):
        event_t=int(path.stem)*.5;target=event_t+offset;row=db.execute('select timestamp,data from messages where topic_id=110 and timestamp>=? order by timestamp limit 1',(round(target*1e9),)).fetchone()
        if not row:continue
        ir,_=parse_image(row[1]);ev=cv2.imread(str(path));canvas=np.zeros((420,1280,3),np.uint8)
        canvas[60:]=np.hstack([cv2.resize(ev,(640,360)),cv2.cvtColor(ir,cv2.COLOR_GRAY2BGR)])
        label=next((p['label'] for p in phases if p['event_start_s']<=event_t<p['event_end_s']),'end')
        cv2.putText(canvas,f'Event {event_t:.2f}s / sampled 2ms', (10,22),0,.6,(255,255,255),1)
        cv2.putText(canvas,f'D455 db {row[0]/1e9:.2f}s / exposure 8.5ms',(650,22),0,.6,(255,255,255),1)
        cv2.putText(canvas,'ROUGH OFFSET ONLY | '+label,(10,47),0,.6,(0,220,255),1);writer.write(canvas)
        if event_t in [1,5,10,17,25,32,39,47]:contact.append(canvas)
    writer.release();cap=cv2.VideoCapture(str(OUT/'paired_coarse_review.mp4'));count=0
    while True:
        ok,im=cap.read()
        if not ok:break
        assert im.shape==(420,1280,3);count+=1
    cap.release();save('video_check.json',dict(frames_decoded=count,fps=2,sampled_event_exposure_ms=2,frame_step_s=.5,dvs_temporal_coverage_fraction=.004,interpolated=False,alignment='coarse only'))
    cv2.imwrite(str(OUT/'paired_evidence.jpg'),np.vstack([cv2.resize(im,(960,315)) for im in contact]))
    print('DONE',audit,flush=True)
    hashes={}
    for path in [RAW,SESSION/'realsense.db3']:
        h=hashlib.sha256()
        with path.open('rb') as f:
            for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
        hashes[str(path)]=dict(bytes=path.stat().st_size,sha256=h.hexdigest());print('HASHED',path.name,flush=True)
    save('source_manifest.json',hashes)

if __name__=='__main__':main()
