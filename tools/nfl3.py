import json, numpy as np
C3=json.load(open('c3.json'))
def yk(c, x, r=0.5):
    cx,cy=np.array(c[0]),np.array(c[1])
    if x<cx[0]: return x/r*1.0001 if x/r<cy[0] else cy[0]  # before start: on the 0.5 boundary
    if x>=cx[-1]: return cy[-1]
    return float(np.interp(x,cx,cy))
def nfl3(ch, x, y, r=0.5):
    out = x>ch['xmax']+1 or x<r*y-1
    xs=min(x,ch['xmax'])
    L=np.array(ch['levels'],float); Y=np.array([yk(c,xs,r) for c in ch['curves']])
    if y>=Y[0]: return L[0]*(Y[0]/y)**2, True
    if y<=Y[-1]: return L[-1]*(Y[-1]/y)**2, True
    k=np.nonzero(Y>=y)[0].max()
    f=(np.log(y)-np.log(Y[k]))/(np.log(Y[k+1])-np.log(Y[k]))
    return float(np.exp(np.log(L[k])+f*(np.log(L[k+1])-np.log(L[k])))), out
