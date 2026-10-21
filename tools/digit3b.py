import sys, json, numpy as np
from PIL import Image
from scipy import ndimage as nd
from axes import mapping
from digit3 import fit, trace
LEV={15:[.5,.75,1,1.25,1.5,2,2.5],16:[.5,.75,1,1.25,1.5,2,2.5],17:[.5,.75,1,1.25,1.5,2,2.5],
 18:[.5,.75,1,1.25,1.5,2,2.5,3],19:[.5,.75,1,1.25,1.5,2,2.5,3],20:[.5,.75,1,1.25,1.5,2,2.5,3],
 21:[.75,1,1.25,1.5,2,2.5,3,4],22:[.75,1,1.25,1.5,2,2.5,3,4],23:[.75,1,1.25,1.5,2,2.5,3,4],
 24:[1.5,2,2.5,3,4,5,7,10],25:[2,2.5,3,4,5,7,10],26:[2.5,3,4,5,7,10]}
def run(n, levels=None, ratio=0.5):
    fn='../ref/nflrepo/A1_%02d.png'%n
    img=Image.open(fn).convert('L'); arr=np.array(img)
    m=mapping(fn); sx,cx,rx=fit(m['bx']); sy,cy,ry=fit(m['ry'])
    left,right,top,bottom=m['left'],m['right'],m['top'],m['bottom']
    bold=nd.binary_opening(arr<100,structure=np.ones((2,2)))
    # remove frame
    bold[:, :left+3]=0; bold[bottom-2:,:]=0; bold[:top+2,:]=0; bold[:, right-2:]=0
    # bold vertical gridlines: columns with long vertical runs
    vgrid=np.zeros(bold.shape[1],bool)
    for c in range(bold.shape[1]):
        col=bold[:,c]; 
        if col.sum()>0:
            lab,nr=nd.label(col); 
            if nr and nd.sum(col,lab,range(1,nr+1)).max()>90: vgrid[c]=True
    # band right of steep diagonal x_mm = ratio*y_mm
    H,W=bold.shape; starts=[]
    yy,xx=np.nonzero(bold)
    xmm=(xx-cx)/sx; ymm=(yy-cy)/sy
    dist=xmm-ratio*ymm
    sel=(dist>0)&(dist<120)&(~vgrid[xx])&(ymm>50)
    pts=np.c_[xx[sel],yy[sel]]
    # cluster by y (connected components in band mask)
    band=np.zeros_like(bold); band[pts[:,1],pts[:,0]]=1
    band=nd.binary_dilation(band,iterations=1)
    lab,nl=nd.label(band)
    for i in range(1,nl+1):
        ys,xs=np.nonzero(lab==i)
        if len(xs)<6: continue
        k=np.argmin(ys); starts.append((int(xs[k]),int(ys[k]),len(xs)))
    starts.sort(key=lambda s:s[1])
    curves=[]
    for x0,y0,_ in starts:
        p=trace(bold,x0,y0,right)
        if len(p)<15: continue
        curves.append(p)
    # dedupe curves ending at same flat y
    ends=[]; uniq=[]
    for p in curves:
        e=(round(p[-1][1]),round(p[-1][0]))
        if any(abs(e[0]-u[0])<3 and abs(e[1]-u[1])<6 for u in ends): continue
        ends.append(e); uniq.append(p)
    uniq.sort(key=lambda p:p[len(p)//2][1] if False else p[-1][1])
    return fn,(sx,cx,sy,cy),uniq,(left,right,top,bottom)
if __name__=='__main__':
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    out={}
    for n in [int(a) for a in sys.argv[1:]]:
        fn,mp,cv,fr=run(n); lev=LEV[n]
        im=Image.open(fn).convert('RGB'); W,H=im.size
        fig=plt.figure(figsize=(W/100,H/100),dpi=100); ax=fig.add_axes([0,0,1,1]); ax.imshow(im); ax.axis('off')
        for i,p in enumerate(cv):
            xs=[q[0] for q in p]; ys=[q[1] for q in p]; ax.plot(xs,ys,'-',color='red' if i%2==0 else 'blue',lw=1)
            ax.text(xs[-1]+2,ys[-1],str(lev[i]) if i<len(lev) else '?',color='red',fontsize=7)
        fig.savefig('dg%02d.png'%n,dpi=100); plt.close(fig)
        print(n,'curves',len(cv),'levels',len(lev))
        out[n]={'map':mp,'levels':lev,'curves':[[((x-mp[1])/mp[0],(y-mp[3])/mp[2]) for x,y in p] for p in cv]}
    json.dump(out,open('d3_%s.json'%'_'.join(sys.argv[1:]),'w'))
LEV.update({36:[.5,.75,1,1.25,1.5,2,2.5,3],37:[.5,.75,1,1.25,1.5,2,2.5,3],38:[.75,1,1.25,1.5,2,2.5,3,4,5],
 39:[.75,1,1.25,1.5,2,2.5,3,4],40:[.75,1,1.25,1.5,2,2.5,3,4],41:[1,1.5,2,2.5,3,4,5,7,10],42:[2,2.5,3,4,5,7,10]})
