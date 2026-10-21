import sys, json, numpy as np
from digit3e import run, show
from starts4 import LV
res={}
for a in sys.argv[1:]:
    n=int(a); fn='../ref/nflrepo/A1_%02d.png'%n; lev=LV[n]
    r=run(fn,len(lev),ratio=1.0,band_mm=(-20,350))
    show(fn,r,lev,'d4l%02d.png'%n)
    print(n,'curves',len(r['curves']),'levels',len(lev),flush=True)
    res[n]={'map':[float(v) for v in r['map']],'levels':lev,'curves':[P.round(1).tolist() for P in r['curves']]}
json.dump(res,open('d4l.json','w'))
