import sys, json, numpy as np
from PIL import Image
from scipy import ndimage as nd
from axes import mapping
from digit3 import fit, trace
from digit3b import LEV
def run(fn, nlev, ratio=0.5, band_mm=(-20,450), th=100):
    arr=np.array(Image.open(fn).convert('L')).astype(float)
    m=mapping(fn); sx,cx,rx=fit(m['bx']); sy,cy,ry=fit(m['ry'])
    left,right,top,bottom=m['left'],m['right'],m['top'],m['bottom']
    bold=nd.binary_opening(arr<th,structure=np.ones((2,2)))
    bold[:, :left+3]=0; bold[bottom-2:,:]=0; bold[:top+2,:]=0; bold[:, right-2:]=0
    H,W=bold.shape
    vg=np.zeros(W,bool)
    for c in range(W):
        lab,nr=nd.label(bold[:,c])
        if nr and nd.sum(bold[:,c],lab,range(1,nr+1)).max()>45: vg[c]=True
    yy,xx=np.nonzero(bold); xmm=(xx-cx)/sx; ymm=(yy-cy)/sy
    d=xmm-ratio*ymm; sel=(d>band_mm[0])&(d<band_mm[1])&(ymm>30)&(~vg[xx])
    seeds=list(zip(xx[sel],yy[sel]))
    seeds=seeds[::2]
    traces=[]
    for x,y in seeds:
        p=np.array(trace(bold,int(x),int(y),right))
        if len(p)<20: continue
        k=0
        while k<len(p) and (vg[int(round(p[k][0]))] or vg[min(W-1,int(round(p[k][0]))+1)] or vg[max(0,int(round(p[k][0]))-1)]): k+=1
        p=p[k:]
        if len(p)<20: continue
        P=np.array([((a-cx)/sx,(b-cy)/sy) for a,b in p])
        if P[-1,0]-P[0,0]<300: continue
        if P[0,1]-P[-1,1] < 0.04*P[-1,1]: continue
        if P[0,0]-ratio*P[0,1] > 300: continue
        traces.append(P)
    traces.sort(key=lambda P:-(P[-1,0]-P[0,0])); uniq=[]
    for P in traces:
        fy=np.median(P[-15:,1])
        if any(abs(fy-np.median(Q[-15:,1]))<max(0.025*fy,8) for Q in uniq): continue
        uniq.append(P)
    uniq.sort(key=lambda P:-np.median(P[-15:,1]))
    return dict(map=(sx,cx,sy,cy),resid=(rx,ry),frame=(left,right,top,bottom),curves=uniq)
def show(fn,r,lev,out):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    im=Image.open(fn).convert('RGB'); W,H=im.size; sx,cx,sy,cy=r['map']
    fig=plt.figure(figsize=(W/100,H/100),dpi=100); ax=fig.add_axes([0,0,1,1]); ax.imshow(im); ax.axis('off')
    for i,P in enumerate(r['curves']):
        ax.plot(P[:,0]*sx+cx,P[:,1]*sy+cy,'-',color=['red','blue'][i%2],lw=1)
        ax.text(P[-1,0]*sx+cx+2,P[-1,1]*sy+cy,str(lev[i]) if i<len(lev) else '?',color='red',fontsize=7)
    fig.savefig(out,dpi=100); plt.close(fig)
if __name__=='__main__':
    res={}
    for a in sys.argv[1:]:
        n=int(a); fn='../ref/nflrepo/A1_%02d.png'%n; lev=LEV[n]
        r=run(fn,len(lev)); show(fn,r,lev,'dg%02d.png'%n)
        print(n,'curves',len(r['curves']),'levels',len(lev),'resid',[round(float(v),1) for v in r['resid']],flush=True)
        res[n]={'map':[float(v) for v in r['map']],'levels':lev,'curves':[P.round(1).tolist() for P in r['curves']]}
    json.dump(res,open('d3e.json' if len(sys.argv)>2 else 'd3e_%s.json'%sys.argv[1],'w'))
