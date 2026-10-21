"""Build GL300_new.html (usage: build.py tools/data GL300_new.html): inject NFL tables into src/index.src.html.
usage: python3 build.py <tables_dir> <out_html>
"""
import json, sys, glob, os, math
import numpy as np

tdir, out = sys.argv[1], sys.argv[2]
src = open(os.path.join(os.path.dirname(__file__), '..', 'src', 'index.src.html'), encoding='utf-8').read()

def comb(J1, J2):
    J1 = np.array(J1, float); J2 = np.array(J2, float)
    return np.logaddexp(J1, J2)

def plate_tab(files):
    rows = []
    for f in files:
        d = json.load(open(f))
        for i, ar in enumerate(d['AR'][:len(d['J'])]):
            rows.append((ar, comb(d['J'][i], d['J2'][i]), np.array(d['w'][i]), np.log(d['qhat'])))
    rows.sort(key=lambda r: r[0])
    lq = rows[0][3]
    return {'AR': [r[0] for r in rows], 'lq': [round(float(x), 5) for x in lq],
            'J': [[round(float(x), 4) for x in r[1]] for r in rows],
            'w': [[float('%.5g' % x) for x in r[2]] for r in rows]}

data = {'t4': plate_tab([os.path.join(tdir, 'tab4.json')]),
        't3': plate_tab([os.path.join(tdir, 'tab3.json')])}
# E1300 chart contours (digitised): 3-sided monolithic/laminated, laminated 4-sided
data['c3'] = json.load(open(os.path.join(tdir, 'c3.json')))
data['c4l'] = json.load(open(os.path.join(tdir, 'c4l.json')))
lam = {}
for f in sorted(glob.glob(os.path.join(tdir, 'lam4_*.json'))):
    d = json.load(open(f))
    rd = lambda v, n: None if v is None else round(v, n)
    lam[d['des']] = {'ply': d['ply'], 'sizes': d['sizes'], 'AR': d['AR'][:len(d['lnB'])], 'q': [round(x, 5) for x in d['q']],
                     'lnB': [[[rd(v, 4) for v in c] for c in row] for row in d['lnB']],
                     'w': [[[None if v is None else float('%.5g' % v) for v in c] for c in row] for row in d['w']]}
data['lam4'] = lam
# Fig. A1.43 laminated two-sided: NFL = C / L^2 (kPa, mm), digitised
data['lam2C'] = {6: 474242, 8: 903649, 10: 1374000, 12: 2360000, 16: 3846000, 19: 5559000}
js = json.dumps(data, separators=(',', ':'))
html = src.replace('/*DATA*/null/*END*/', js)
assert html != src
open(out, 'w', encoding='utf-8').write(html)
print('ok', len(html) // 1024, 'KB', 'lam designations:', sorted(lam))
