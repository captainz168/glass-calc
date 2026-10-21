import json, numpy as np
from axes import mapping
from digit3 import fit
from starts4 import LV
d=json.load(open('d4l.json')); out={}
for n,v in d.items():
    lev=LV[int(n)]; C=[np.array(P) for P in v['curves']]
    m=mapping('../ref/nflrepo/A1_%02d.png'%int(n)); sy,cy,_=fit(m['ry']); ymax=(m['top']-cy)/sy
    good=[]
    for P in C:
        x0,y0=P[0]
        if abs(x0-y0)<0.03*x0+15 or abs(x0-600)<25 or abs(y0-ymax)<40:
            if any(np.hypot(*(P[0]-Q[0]))<0.03*(x0+y0)/2+10 for Q in good): continue
            good.append(P)
    good.sort(key=lambda P:-((P[0][0]+P[0][1]) if abs(P[0][1]-ymax)>=40 else 1e5+P[0][0]))
    print(n,'valid',len(good),'levels',len(lev),[ (round(P[0][0]),round(P[0][1])) for P in good])
    out[n]={'map':v['map'],'levels':lev[:len(good)],'curves':[P.round(1).tolist() for P in good]}
json.dump(out,open('d4lclean.json','w'))
