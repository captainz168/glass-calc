import json, itertools, numpy as np, sys
d=json.load(open(sys.argv[1])); out={}
for n,v in d.items():
    lev=v['levels']; C=[np.array(P) for P in v['curves']]
    fy=np.array([np.median(P[-15:,1]) for P in C])
    best=None
    for comb in itertools.combinations(range(len(C)),len(lev)):
        vals=np.log(np.array(lev)*fy[list(comb)]**2)
        s=vals.std()
        if best is None or s<best[0]: best=(s,comb)
    s,comb=best
    out[n]={'map':v['map'],'levels':lev,'curves':[C[i].round(1).tolist() for i in comb]}
    print(n,'std %.3f'%s,'dropped',[int(fy[i]) for i in range(len(C)) if i not in comb],'products',[round(float(lev[k]*fy[i]**2/1e6),3) for k,i in enumerate(comb)])
json.dump(out,open(sys.argv[2],'w'))
