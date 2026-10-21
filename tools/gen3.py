import numpy as np, json, sys
from jtab import jcurve
part=int(sys.argv[1])
qh=np.geomspace(0.05,1e4,40)
ARs=[0.5,0.6,0.75,0.9,1.0,1.25,1.5,1.75,2.0,2.5,3.0,4.0,5.0,6.0,8.0,10.0,14.0,20.0]
mine=ARs[part::2]
out={'qhat':qh.tolist(),'AR':mine,'J':[],'J2':[],'w':[]}
for AR in mine:
    r=jcurve(AR,'LRB',qh,nw=18,nuv=16,ng=28)
    out['J'].append(r[:,1].tolist()); out['J2'].append(r[:,2].tolist()); out['w'].append(r[:,3].tolist())
    json.dump(out,open('tab3_%d.json'%part,'w')); print(AR,flush=True)
