import numpy as np, json, sys
from jtab import jcurve
qh=np.geomspace(0.02,6e4,70)
ARs=[1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.8,2.0,2.25,2.5,2.75,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0]
out={'qhat':qh.tolist(),'AR':ARs,'J':[],'J2':[],'w':[]}
for AR in ARs:
    r=jcurve(AR,'LRBT',qh,nw=28,nuv=26,ng=40)
    out['J'].append(r[:,1].tolist()); out['J2'].append(r[:,2].tolist()); out['w'].append(r[:,3].tolist())
    print(AR,flush=True)
json.dump(out,open('tab4.json','w'))
