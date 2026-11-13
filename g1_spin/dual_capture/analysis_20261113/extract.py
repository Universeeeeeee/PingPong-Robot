import sys,json,sqlite3
from pathlib import Path
import numpy as np,cv2
from audit_remote import OUT,SESSION,RAW,parse_image,save
cv2.setNumThreads(1)

def color():
 db=sqlite3.connect('file:'+str(SESSION/'color.db3')+'?mode=ro',uri=True);rough=json.loads((OUT/'rough_alignment.json').read_text())['corresponding_db_s'];rows=[];crops=[]
 for ts,b in db.execute('select timestamp,data from messages where topic_id=102 and timestamp>=? order by timestamp',(int((rough-4)*1e9),)):
  im,h=parse_image(b);hsv=cv2.cvtColor(im,cv2.COLOR_BGR2HSV);mask=cv2.inRange(hsv,(0,110,65),(30,255,255));mask[:90]=0;mask[260:]=0;mask[:,:180]=0;mask[:,410:]=0
  cs,_=cv2.findContours(mask,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE);cs=[c for c in cs if 250<cv2.contourArea(c)<1800]
  if not cs:continue
  c=max(cs,key=cv2.contourArea);(x,y),r=cv2.minEnclosingCircle(c)
  if not 10<r<28:continue
  # Same normalized ball image at every frame; preserve color contrast for periodicity.
  M=np.float32([[r/26,0,x-r*32/26],[0,r/26,y-r*32/26]])
  crop=cv2.warpAffine(im,M,(64,64),flags=cv2.INTER_LINEAR|cv2.WARP_INVERSE_MAP)
  rows.append([ts/1e9-rough,x,y,r]);crops.append(crop)
 np.savez_compressed(OUT/'color_track.npz',track=np.array(rows),crops=np.array(crops));print('COLOR',len(rows),flush=True)

def events():
 sys.path.insert(0,str(Path.home()/'g1_spin/scripts'));from event_io import source
 info,batches=source(RAW,batch_size=2000000);origin=info['sdk_index_start_us'];track=[];imgs={};windows={k:[] for k in [9,17,19,25,31,33,42,45]}
 for e in batches:
  t=(e['t'].astype(np.int64)-origin);roi=(e['x']>200)&(e['x']<850)&(e['y']>100)&(e['y']<550);e=e[roi];t=t[roi]
  for k in windows:
   z=e[(t>=k*1000000)&(t<(k+.25)*1000000)]
   if len(z):windows[k].append(z)
  select=t%100000<20000;ss=t[select]//100000;ee=e[select]
  for k in np.unique(ss):
   z=ee[ss==k];im=imgs.setdefault(int(k),np.zeros((450,650),np.uint8));im[z['y']-100,z['x']-200]=255
 for k,im in sorted(imgs.items()):
  cc=cv2.HoughCircles(cv2.GaussianBlur(im,(5,5),1),cv2.HOUGH_GRADIENT,1,80,param1=60,param2=25,minRadius=40,maxRadius=85)
  candidates=[]
  if cc is not None:
   for x,y,r in cc[0]:
    if y+r>440:continue
    candidates.append([float(x+200),float(y+100),float(r)])
  track.append(dict(t=k*.1+.01,candidates=candidates))
 save('event_track.json',track);(OUT/'windows').mkdir(exist_ok=True)
 for k,parts in windows.items():np.savez_compressed(OUT/f'windows/t{k}.npz',events=np.concatenate(parts),origin=origin)
 print('EVENT_TRACK',sum(len(r['candidates'])==1 for r in track),len(track),flush=True)
if __name__=='__main__':color() if sys.argv[1]=='color' else events()
