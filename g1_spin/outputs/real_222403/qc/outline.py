exec(open('/Users/vae/Projects/PingPongRobot/g1_spin/outputs/real_222403/qc/survey.py').read().split('counts=np.bincount')[0])
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.optimize import least_squares
rows=[]
fig,axs=plt.subplots(3,4,figsize=(14,10))
selected=[.1,1,2,3,4,6,8,10,12,14,15,16]
for t in np.arange(.1,16.2,.1):
 a,b=np.searchsorted(ts,[round(t*1e6),round((t+.1)*1e6)])
 s=e[a:b];c=np.bincount(s['y'].astype(int)*w+s['x'],minlength=h*w).reshape(h,w)
 smooth=gaussian_filter(c.astype(float),1)
 angles=np.linspace(-np.pi,0,121); radii=np.arange(30,76,.5)
 xx=556+np.cos(angles[:,None])*radii;yy=224+np.sin(angles[:,None])*radii
 profiles=map_coordinates(smooth,[yy,xx],order=1)
 fits=[];points=None
 for threshold in [.3,.6,1.]:
  rr=[];theta=[]
  for angle,p in zip(angles,profiles):
   indices=np.where(p>=threshold)[0]
   if len(indices):rr.append(radii[indices[-1]]);theta.append(angle)
  if len(rr)<45:continue
  points=np.column_stack((556+np.cos(theta)*rr,224+np.sin(theta)*rr))
  def residual(p):return np.linalg.norm(points-p[:2],axis=1)-p[2]
  fit=least_squares(residual,[556,224,50],loss='soft_l1',f_scale=1,bounds=([525,200,35],[590,250,70]))
  rms=float(np.sqrt(np.mean(residual(fit.x)**2)))
  fits.append([*fit.x,rms,len(points)])
 if not fits:continue
 best=fits[min(1,len(fits)-1)]
 rows.append([t,*best,*np.std(np.array(fits)[:,:3],axis=0)])
 if any(abs(t-v)<.001 for v in selected):
  idx=min(range(len(selected)),key=lambda i:abs(selected[i]-t));ax=axs.flat[idx]
  ax.imshow(np.log1p(c),cmap='magma',extent=[0,w,h,0],vmin=0,vmax=3)
  cx,cy,r=best[:3];ax.add_patch(plt.Circle((cx,cy),r,fill=False,color='cyan',lw=1));ax.plot(cx,cy,'c+')
  ax.set_xlim(475,640);ax.set_ylim(280,150);ax.set_title(f'{t:.1f}s: {cx:.1f},{cy:.1f},r={r:.1f}; rms={best[3]:.1f}')
fig.tight_layout();fig.savefig(out/'circle_fits.png',dpi=150);plt.close(fig)
a=np.array(rows)
np.savetxt(out/'circle_fits.csv',a,delimiter=',',header='relative_s,cx,cy,radius,rms,arc_points,cx_threshold_std,cy_threshold_std,r_threshold_std',comments='')
fig,axs=plt.subplots(3,1,figsize=(12,7),sharex=True)
for ax,col,title in zip(axs,[1,2,3],['Circle cx','Circle cy','Radius']):ax.plot(a[:,0],a[:,col],'.-');ax.set_ylabel(title);ax.grid(alpha=.3)
axs[-1].set_xlabel('Relative seconds');fig.tight_layout();fig.savefig(out/'circle_drift.png',dpi=150)
for start,end in [(0,3),(3,6),(6,9),(9,12),(12,14),(14,16)]:
 sel=a[(a[:,0]>=start)&(a[:,0]<end)&(a[:,4]<4)]
 print(start,end,'n',len(sel),'median',np.median(sel[:,1:4],axis=0),'range',np.ptp(sel[:,1:4],axis=0),'median_rms',np.median(sel[:,4]))
