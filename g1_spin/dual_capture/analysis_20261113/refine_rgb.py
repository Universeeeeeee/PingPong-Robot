import json
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
P=Path(__file__).resolve().parent
rows=json.loads((P/'rgb_periodicity.json').read_text());out=[]
for row in rows:
 t=np.array(row['lags']);s=np.array(row['scores']);peaks,_=find_peaks(s,prominence=.15);fit=[]
 for i in peaks:
  if not 0<i<len(t)-1:continue
  poly=np.polyfit(t[i-1:i+2]-t[i],s[i-1:i+2],2);tp=t[i]-poly[1]/(2*poly[0])
  if .20<tp<1.08:fit.append((tp,s[i]))
 if not fit:continue
 base=min(x[0] for x in fit);kept=[(round(tp/base),tp,c) for tp,c in fit if abs(tp/base-round(tp/base))<.15];n=np.array([x[0] for x in kept]);v=np.array([x[1] for x in kept]);period=float(n@v/(n@n))
 out.append(dict(begin=row['begin'],end=row['end'],rpm=60/period,period_s=period,individual_rpm=[60*k/t for k,t,c in kept],peaks=kept,note='Subframe autocorrelation peak interpolation and integer-multiple period fit. Not a calibrated accuracy interval.'))
(P/'rgb_period_refined.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
