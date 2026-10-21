import sys, json, numpy as np
from PIL import Image
from scipy import ndimage as nd
from axes import mapping
from digit3 import fit
LV={29:[.5,.75,1,1.25,1.5,2,2.5,3,4,5,7,10],30:[.5,.75,1,1.25,1.5,2,2.5,3,4,5,7,10],31:[.75,1,1.25,1.5,2,2.5,3,4,5,7,10],
    32:[1,1.25,1.5,2,2.5,3,4,5,7,10],33:[1.25,1.5,2,2.5,3,4,5,7,10],34:[1.5,2,2.5,3,4,5,7,10],35:[2,2.5,3,4,5,7,10],
    7:[.5,.75,1,1.25,1.5,2,2.5,3,4,5,7,10,15]}
def starts(n):
    fn='../ref/nflrepo/A1_%02d.png'%n; arr=np.array(Image.open(fn).convert('L')).astype(float)
    m=mapping(fn); sx,cx,_=fit(m['bx']); sy,cy,_=fit(m['ry'])
    bold=nd.binary_opening(arr<100,structure=np.ones((2,2)))
    # walk along AR=1 line, probe a point offset perpendicular (down-right) by 4 px
    ss=np.arange(300,4000,2.0); hits=[]
    nx,ny=1/np.sqrt(2),1/np.sqrt(2)  # perpendicular in px (down-right) approx
    for s in ss:
        px=cx+sx*s; py=cy+sy*s
        found=False
        for off in (3,4,5,6):
            qx,qy=int(round(px+off*0.7)),int(round(py+off*0.7))
            if 0<=qy<bold.shape[0] and 0<=qx<bold.shape[1] and bold[qy,qx]: found=True
        hits.append(found)
    hits=np.array(hits); lab,k=nd.label(hits)
    cen=[float(ss[lab==i].mean()) for i in range(1,k+1)]
    return cen, LV[n]
if __name__=='__main__':
    n=int(sys.argv[1]); c,l=starts(n); print(n,'starts',[round(x) for x in c],'levels',l[::-1])
