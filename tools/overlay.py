"""Overlay computed NFL contours on ASTM chart image.
usage: overlay.py chart.png tab.json tmin x0px xscale(px/mm) y0px yscale(px/mm) levels out.png [mode]
mode '4': x=long, y=short (a=short,b=long) ; mode '3': x=parallel supported edge, y=free edge
"""
import sys, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from nflcalc import Table
fn, tab, tmin, x0, xs, y0, ys, levels, out = sys.argv[1:10]
mode = sys.argv[10] if len(sys.argv) > 10 else '4'
tmin = float(tmin); x0, xs, y0, ys = map(float, (x0, xs, y0, ys))
levels = [float(v) for v in levels.split(',')]
T = Table(tab)
im = Image.open(fn).convert('RGB'); W, H = im.size
X = np.linspace(100, 6500, 240); Y = np.linspace(100, 4000, 200)
Z = np.full((len(Y), len(X)), np.nan)
for j, y in enumerate(Y):
    for i, x in enumerate(X):
        if mode == '4':
            if y > x: continue
            Z[j, i] = T.nfl(y, x, tmin)
        else:
            if x < 0.5 * y: continue
            Z[j, i] = T.nfl(y, x, tmin, AR=x / y)  # free edge = y (x-length), supported = x
fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(im); ax.axis('off')
PX = x0 + xs * X; PY = y0 - ys * Y
cs = ax.contour(PX, PY, Z, levels=levels, colors='red', linewidths=1.0)
ax.clabel(cs, fontsize=6, fmt='%.2f', colors='red')
ax.set_xlim(0, W); ax.set_ylim(H, 0)
fig.savefig(out, dpi=100)
