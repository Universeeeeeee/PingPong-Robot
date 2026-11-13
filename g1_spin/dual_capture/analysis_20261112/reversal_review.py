import sqlite3,json
import cv2,numpy as np
from audit_remote import OUT,SESSION,parse_image,save
rough=json.loads((OUT/'rough_alignment.json').read_text())['corresponding_db_s'];db=sqlite3.connect('file:'+str(SESSION/'realsense.db3')+'?mode=ro',uri=True)
frames=[];records=[]
for event_t in np.arange(25,36,.25):
 target=event_t+rough;row=db.execute('select timestamp,data from messages where topic_id=110 and timestamp>=? order by timestamp limit 1',(round(target*1e9),)).fetchone();im,_=parse_image(row[1])
 # Enlarge the ball, drill and hands; does not add temporal samples.
 crop=im[70:350,190:500];tile=cv2.resize(crop,(310,280));canvas=np.zeros((310,310),np.uint8);canvas[30:]=tile;cv2.putText(canvas,f'~event {event_t:.2f}s',(5,21),0,.5,255,1);frames.append(canvas)
 records.append(dict(event_coarse_s=float(event_t),db_s=row[0]/1e9))
for page in range((len(frames)+11)//12):
 canvas=np.zeros((930,1240),np.uint8)
 for j,frame in enumerate(frames[page*12:(page+1)*12]):
  y,x=divmod(j,4);canvas[y*310:(y+1)*310,x*310:(x+1)*310]=frame
 cv2.imwrite(str(OUT/f'reversal_detail_{page}.jpg'),canvas)
save('reversal_detail_times.json',records)
