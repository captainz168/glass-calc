import numpy as np, json, sys, time
from lamplate import LamPlate
from nflcalc import K, DUR, TARGET
DES={5:2.16,6:2.92,8:3.78,10:4.57,12:5.56,16:7.42,19:9.02}
des=[int(x) for x in sys.argv[1].split(',')]
sizes=np.geomspace(250,3500,12).round().tolist(); ARs=[1.0,1.25,1.6,2.0,3.0,5.0]
qs=np.geomspace(0.1,40,22)
for d in des:
    h=DES[d]; out={'des':d,'ply':h,'pvb':0.76,'G':0.44e6,'sizes':sizes,'AR':ARs,'q':qs.tolist(),'lnB':[],'w':[]}
    for AR in ARs:
        rB=[];rW=[]
        for a in sizes:
            p=LamPlate(a*1e-3,a*AR*1e-3,h*1e-3,h*1e-3,0.76e-3,0.44e6,nw=12,nuv=10,ng=20)
            lB=[];W=[];qn=None
            for q in qs:
                if qn is None and lB and lB[-1]>np.log(TARGET): qn=q
                if qn is not None and q>5.5*qn: lB.append(None); W.append(None); continue
                p.solve(q*1e3); lB.append(float(np.log(K*DUR*sum(p.risk())))); W.append(float(p.wmax()*1e3))
            rB.append(lB); rW.append(W)
        out['lnB'].append(rB); out['w'].append(rW); print(d,AR,flush=True)
        json.dump(out,open('lam4_%d.json'%d,'w'))
