import numpy as np, sys, json, time
from plate import *
m=7
def jcurve(AR, edges, qhats, nw=16, nuv=14, ng=26, ab=1e4, along='a'):
    """AR = b/a with a = x-length. returns J(qhat) per face and w/t max"""
    a=np.sqrt(ab/AR); b=a*AR
    p=Plate(a,b,edges,nw=nw,nuv=nuv,ng=ng)
    out=[]
    for qh in qhats:
        p.solve(qh/ab**2)
        S=risk_integral(p)
        J=[np.log(s*ab**m/ab) if s>0 else -np.inf for s in S]
        out.append((qh,J[0],J[1],p.wmax()))
    return np.array(out)
if __name__=='__main__':
    qh=np.geomspace(0.1,3e5,8)
    for AR in [1.0,5.0]:
        for nw,nuv,ng in [(16,14,26),(22,20,34),(30,28,44)]:
            t0=time.time(); r=jcurve(AR,'LRBT',qh,nw,nuv,ng)
            print(AR,nw,' '.join('%.3f'%x for x in r[:,1]),'| w/t',' '.join('%.2f'%x for x in r[:,3]),'%.1fs'%(time.time()-t0)); sys.stdout.flush()
