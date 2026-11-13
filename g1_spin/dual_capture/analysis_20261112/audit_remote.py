import argparse,json,sqlite3,struct,sys,time,hashlib
from pathlib import Path
import numpy as np
import cv2
HOME=Path.home();SESSION=HOME/'g1_spin/dual_capture/recordings/formal_20261112_01';RAW=HOME/'g1_spin/data/raw/DSP03F100025_20261112_224226.raw';OUT=HOME/'g1_spin/dual_capture/analysis_20261112';OUT.mkdir(exist_ok=True)

def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def parse_image(b):
    if b[:2]!=b'\x00\x01':raise ValueError('Only CDR little endian supported')
    offset=4
    def align(n):
        nonlocal offset
        offset=4+((offset-4+n-1)//n)*n
    def u32():
        nonlocal offset
        align(4);x=struct.unpack_from('<I',b,offset)[0];offset+=4;return x
    sec=u32();ns=u32();length=u32();frame=b[offset:offset+length-1].decode();offset+=length
    h=u32();w=u32();length=u32();encoding=b[offset:offset+length-1].decode();offset+=length
    endian=b[offset];offset+=1;step=u32();length=u32();payload=b[offset:offset+length]
    assert len(payload)==length and encoding=='8UC1' and endian==0 and length==h*step
    return np.frombuffer(payload,np.uint8).reshape(h,step)[:,:w].copy(),dict(header_s=sec+ns/1e9,frame=frame,width=w,height=h)

def realsense():
    db=sqlite3.connect('file:'+str(SESSION/'realsense.db3')+'?mode=ro',uri=True)
    topics=db.execute('select id,name,type from topics').fetchall();counts=[]
    # One indexed pass over messages for counts and timestamp extent.
    for tid,count,first,last in db.execute('select topic_id,count(*),min(timestamp),max(timestamp) from messages group by topic_id'):
        name=next(x[1] for x in topics if x[0]==tid)
        if '/image/data' in name:counts.append(dict(topic=name,id=tid,count=count,first_db_ns=first,last_db_ns=last,span_s=(last-first)/1e9))
    save('db_inventory.json',counts)
    mid=next(t[0] for t in topics if 'Infrared_1/image/metadata' in t[1]);iid=next(t[0] for t in topics if 'Infrared_1/image/data' in t[1])
    metadata=[]
    for dbts,b in db.execute('select timestamp,data from messages where topic_id=? order by timestamp',(mid,)):
        length=struct.unpack_from('<I',b,4)[0];text=b[8:8+length-1].decode();items=dict(x.split('=',1) for x in text.split(';') if '=' in x)
        metadata.append(dict(db_ns=dbts,**items))
    save('ir_left_metadata.json',metadata)
    nums=np.array([int(x['Frame number']) for x in metadata]);gaps=np.diff(nums)-1
    sensor=np.array([int(x['Sensor Timestamp']) for x in metadata]);domain=[x['timestamp_domain'] for x in metadata]
    save('db_frame_audit.json',dict(frames=len(nums),frame_number_gaps=int(gaps[gaps>0].sum()),duplicate_or_reversed=int((np.diff(nums)<=0).sum()),sensor_timestamp_regressions=int((np.diff(sensor)<0).sum()),sensor_span_s=float((sensor[-1]-sensor[0])*1e-6),sensor_rate_hz=float((len(sensor)-1)*1e6/(sensor[-1]-sensor[0])),domains=list(dict.fromkeys(domain)),exposures=sorted(set(x.get('Actual Exposure') for x in metadata)),callback_count=json.loads((SESSION/'session.json').read_text())['counts']))
    links=json.loads((SESSION/'insight_events.jsonl').read_text().splitlines()[0]);nearest=None
    for line in (SESSION/'frames.jsonl').open():
        x=json.loads(line)
        if x['stream']=='stream.infrared_1':
            if nearest is None or abs(x['host_monotonic_ns']-links['host_monotonic_ns'])<abs(nearest['host_monotonic_ns']-links['host_monotonic_ns']):nearest=x
    matched=next(x for x in metadata if int(x['Frame number'])==nearest['frame_number']);rough=matched['db_ns']/1e9
    save('rough_alignment.json',dict(raw_file_first_seen_host_ns=links['host_monotonic_ns'],nearest_callback=nearest,corresponding_db_s=rough,note='File detection time is a rough anchor, not optical synchronization. RAW start and write buffering are unmeasured.'))
    image_dir=OUT/'ir_samples';image_dir.mkdir(exist_ok=True);samples=[]
    # Entire active period plus 5-second margins, one image per second.
    for target in np.arange(max(0,rough-5),min(counts[0]['last_db_ns']/1e9,rough+62),1):
        hit=db.execute('select timestamp,data from messages where topic_id=? and timestamp>=? order by timestamp limit 1',(iid,round(target*1e9))).fetchone()
        if not hit:continue
        im,header=parse_image(hit[1]);filename=f'{target:.2f}.png';cv2.imwrite(str(image_dir/filename),im)
        samples.append(dict(db_s=hit[0]/1e9,file='ir_samples/'+filename,**header,mean=float(im.mean()),p90=float(np.percentile(im,90))))
    save('ir_samples.json',samples)
    for page in range((len(samples)+11)//12):
        canvas=np.zeros((3*210,4*320),np.uint8)
        for j,item in enumerate(samples[page*12:(page+1)*12]):
            im=cv2.imread(str(OUT/item['file']),0);tile=cv2.resize(im,(320,180));r,c=divmod(j,4);canvas[r*210+30:r*210+210,c*320:(c+1)*320]=tile
            cv2.putText(canvas,f'db {item["db_s"]:.2f}s',(c*320+5,r*210+21),0,.55,255,1)
        cv2.imwrite(str(OUT/f'ir_contact_{page}.jpg'),canvas)
    print('RS_DONE',len(metadata),'frames; rough',rough,flush=True)

def events():
    sys.path.insert(0,str(HOME/'g1_spin/scripts'));from event_io import source,EventStats
    info,batches=source(RAW,batch_size=2000000);save('event_info.json',info)
    stats=EventStats(info['width'],info['height']);origin=info['sdk_index_start_us'];bins=np.zeros(7000,np.int64);snapshots={};anomalies=[];last=None;hot=np.zeros((720,1280),np.int64)
    # Full audit, previews at 0.5-second spacing with 2ms exposure.
    for bi,e in enumerate(batches):
        stats.add(e);t=e['t'].astype(np.int64);rel=t-origin
        rev=np.flatnonzero(np.diff(t)<0)
        if last is not None and t[0]<last:anomalies.append(dict(previous=last,current=int(t[0]),batch=bi))
        anomalies.extend(dict(previous=int(t[i]),current=int(t[i+1]),batch=bi) for i in rev[:20]);last=int(t[-1])
        b=rel//10000;valid=(b>=0)&(b<len(bins));bins+=np.bincount(b[valid],minlength=len(bins))
        hot+=np.bincount(e['y'].astype(int)*1280+e['x'].astype(int),minlength=720*1280).reshape(720,1280)
        slot=rel//500000;mask=(rel%500000<2000)&(rel>=0)&(slot<125)
        for k in np.unique(slot[mask]):
            q=e[mask&(slot==k)];im=snapshots.setdefault(int(k),np.zeros((720,1280,3),np.uint8))
            for pol,color in [(0,(255,170,70)),(1,(65,175,255))]:
                z=q[q['p']==pol];im[z['y'],z['x']]=color
        if bi%40==0:print('EVENT_PROGRESS',bi,stats.n,flush=True)
    save('event_stats.json',stats.report());save('event_timestamp_anomalies.json',anomalies);np.save(OUT/'event_rate_10ms.npy',bins);np.save(OUT/'event_hotmap.npy',hot)
    (OUT/'event_samples').mkdir(exist_ok=True)
    for k,im in sorted(snapshots.items()):cv2.imwrite(str(OUT/f'event_samples/{k:03d}.png'),im)
    keys=sorted(snapshots)
    for page in range((len(keys)+11)//12):
        canvas=np.zeros((3*210,4*320,3),np.uint8)
        for j,k in enumerate(keys[page*12:(page+1)*12]):
            row,col=divmod(j,4);canvas[row*210+30:row*210+210,col*320:(col+1)*320]=cv2.resize(snapshots[k],(320,180));cv2.putText(canvas,f'event {k*.5:.2f}s / 2ms',(col*320+4,row*210+21),0,.5,(255,255,255),1)
        cv2.imwrite(str(OUT/f'event_contact_{page}.jpg'),canvas)
    print('EVENT_DONE',stats.report(),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['rs','events']);args=p.parse_args();realsense() if args.mode=='rs' else events()
