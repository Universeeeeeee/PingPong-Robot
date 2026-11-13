import argparse,json,sqlite3,struct,sys,time,hashlib
from pathlib import Path
import numpy as np
import cv2
HOME=Path.home();SESSION=HOME/'g1_spin/dual_capture/recordings/formal_rgb90_20261113_02';RAW=HOME/'g1_spin/data/raw/DSP03F100025_20261113_192618.raw';OUT=HOME/'g1_spin/dual_capture/analysis_20261113';OUT.mkdir(exist_ok=True)

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
    assert len(payload)==length and encoding in ('8UC1','rgb8','8UC3') and endian==0 and length==h*step
    return (np.frombuffer(payload,np.uint8).reshape(h,step)[:,:w].copy() if encoding=='8UC1' else cv2.cvtColor(np.frombuffer(payload,np.uint8).reshape(h,step)[:,:w*3].reshape(h,w,3),cv2.COLOR_RGB2BGR)),dict(header_s=sec+ns/1e9,frame=frame,width=w,height=h)

def realsense():
    db=sqlite3.connect('file:'+str(SESSION/'color.db3')+'?mode=ro',uri=True)
    topics=db.execute('select id,name,type from topics').fetchall()
    iid=next(t[0] for t in topics if 'Color_0/image/data' in t[1])
    mid=next(t[0] for t in topics if 'Color_0/image/metadata' in t[1])
    meta=[]
    for ts,b in db.execute('select timestamp,data from messages where topic_id=? order by timestamp',(mid,)):
        length=struct.unpack_from('<I',b,4)[0];text=b[8:8+length-1].decode();items=dict(x.split('=',1) for x in text.split(';') if '=' in x)
        meta.append(dict(db_ns=ts,**items))
    save('color_metadata.json',meta)
    session=json.loads((SESSION/'session.json').read_text());frames=[json.loads(x) for x in (SESSION/'frames.jsonl').open()]
    near=min(frames,key=lambda x:abs(x['host_ns']-session['seen_host_ns']))
    match=next(x for x in meta if int(x['Frame number'])==near['frame_number']);rough=match['db_ns']/1e9
    save('rough_alignment.json',dict(corresponding_db_s=rough,warning='File-detection anchor, not calibrated synchronization'))
    (OUT/'color_samples').mkdir(exist_ok=True);samples=[]
    for target in np.arange(rough-4,rough+59,1):
        row=db.execute('select timestamp,data from messages where topic_id=? and timestamp>=? order by timestamp limit 1',(iid,round(target*1e9))).fetchone()
        if not row:continue
        im,h=parse_image(row[1]);fn=f'color_samples/{target:.2f}.png';cv2.imwrite(str(OUT/fn),im);samples.append(dict(db_s=row[0]/1e9,file=fn,**h))
    save('color_samples.json',samples)
    for page in range((len(samples)+11)//12):
        canvas=np.zeros((630,1280,3),np.uint8)
        for j,item in enumerate(samples[page*12:(page+1)*12]):
            r,c=divmod(j,4);im=cv2.imread(str(OUT/item['file']));canvas[r*210+30:r*210+210,c*320:(c+1)*320]=cv2.resize(im,(320,180));cv2.putText(canvas,f"relative {item['db_s']-rough:.2f}s",(c*320+5,r*210+21),0,.55,(255,255,255),1)
        cv2.imwrite(str(OUT/f'color_contact_{page}.jpg'),canvas)
    print('RS_DONE',len(meta),'rough',rough,flush=True)

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
