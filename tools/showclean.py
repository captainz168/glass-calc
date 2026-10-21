import json, sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from PIL import Image
d=json.load(open(sys.argv[1])); pre=sys.argv[2] if len(sys.argv)>2 else 'cl'
for n,v in d.items():
    fn='../ref/nflrepo/A1_%02d.png'%int(n); sx,cx,sy,cy=v['map']
    im=Image.open(fn).convert('RGB'); W,H=im.size
    fig=plt.figure(figsize=(W/100,H/100),dpi=100); ax=fig.add_axes([0,0,1,1]); ax.imshow(im); ax.axis('off')
    for i,P in enumerate(v['curves']):
        P=np.array(P); ax.plot(P[:,0]*sx+cx,P[:,1]*sy+cy,'-',color=['red','blue'][i%2],lw=1.2)
        ax.text(P[-1,0]*sx+cx+3,P[-1,1]*sy+cy+3,str(v['levels'][i]),color='red',fontsize=8,weight='bold')
    fig.savefig('%s%s.png'%(pre,n),dpi=100); plt.close(fig)
