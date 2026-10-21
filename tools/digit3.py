"""Digitise E1300 three-sided NFL charts (contours) by tracing bold curves from their labels.
usage: digit3.py chart.png out.json [levels]
"""
import sys, json, re, subprocess, numpy as np
from PIL import Image
from scipy import ndimage as nd
from axes import mapping

SP = '/tmp/claude-0/-home-claude/fa34c777-b303-5da7-acc5-8854c649bf21/scratchpad/e1300/'

def fit(pairs):
    v = np.array([p[0] for p in pairs], float); px = np.array([p[1] for p in pairs], float)
    A = np.c_[v, np.ones_like(v)]; s, c = np.linalg.lstsq(A, px, rcond=None)[0]
    resid = np.abs(A @ [s, c] - px).max()
    return s, c, resid

def labels(img, box):
    """OCR contour labels like 0.50, 1.25, 10.0 inside plot (left region)."""
    l, t, r, b = box
    crop = img.crop((l, t, r, b)); sc = 4
    crop = crop.resize((crop.width * sc, crop.height * sc), Image.LANCZOS)
    crop.save(SP + '_l.png')
    out = subprocess.run(['tesseract', SP + '_l.png', '-', '--psm', '11', 'tsv'], capture_output=True, text=True).stdout
    res = []
    for line in out.splitlines()[1:]:
        p = line.split('\t')
        if len(p) < 12: continue
        txt = p[11].strip()
        if re.fullmatch(r'\d{1,2}\.\d{1,2}', txt):
            L, T, W, H = map(int, p[6:10])
            res.append((float(txt), l + (L + W) / sc, t + (T + H / 2) / sc, l + L / sc))
    return res

def trace(bold, x0, y0, xmax, step=2.0):
    """follow bold curve from (x0,y0) heading down, then right."""
    H, W = bold.shape
    x0, y0 = float(x0), float(y0); pts = [(x0, y0)]; ang = np.pi / 2  # pointing down (image y)
    x, y = float(x0), float(y0)
    for it in range(4000):
        best = None
        for d in np.linspace(-0.7, 0.7, 29):
            a = ang + d
            # forbid going up or left
            if np.sin(a) < -0.05 or np.cos(a) < -0.05: continue
            cx, cy = x + step * np.cos(a), y + step * np.sin(a)
            ix, iy = int(round(cx)), int(round(cy))
            if not (1 <= ix < W - 1 and 1 <= iy < H - 1): continue
            sc = bold[iy - 1:iy + 2, ix - 1:ix + 2].sum()
            score = sc - 3.0 * abs(d)
            if sc > 0 and (best is None or score > best[0]): best = (score, a, cx, cy)
        if best is None: break
        _, a, x, y = best
        ang = 0.6 * ang + 0.4 * a
        pts.append((x, y))
        if x >= xmax - 1: break
    return pts

def digitise(fn, levels=None):
    img = Image.open(fn).convert('L'); arr = np.array(img)
    m = mapping(fn)
    sx, cx, rx = fit(m['bx']); sy, cy, ry = fit(m['ry'])
    bold = nd.binary_opening(arr < 100, structure=np.ones((2, 2)))
    left, right, top, bottom = m['left'], m['right'], m['top'], m['bottom']
    labs = labels(img, (left + 2, top + 2, left + int((right - left) * 0.45), bottom - 5))
    labs.sort(key=lambda r: r[0])
    if levels: labs = [l for l in labs if l[0] in levels]
    curves = []
    for (L, lx, ly, lleft) in labs:
        # nearest bold pixel to the right of the label within a window
        ys, xs = np.nonzero(bold[int(ly) - 12:int(ly) + 14, int(lx):int(lx) + 30])
        if len(xs) == 0: continue
        k = np.argmin((xs) ** 2 + (ys - 12) ** 2 * 0.5)
        x0, y0 = int(lx) + xs[k], int(ly) - 12 + ys[k]
        pts = trace(bold, x0, y0, right)
        mm = [((px - cx) / sx, (py - cy) / sy) for px, py in pts]
        curves.append({'L': L, 'px': pts, 'mm': mm})
    return {'file': fn, 'map': [float(sx), float(cx), float(sy), float(cy)], 'resid': [float(rx), float(ry)], 'frame': [left, right, top, bottom],
            'xmax_mm': (right - cx) / sx, 'curves': curves}

if __name__ == '__main__':
    d = digitise(sys.argv[1])
    json.dump(d, open(sys.argv[2], 'w'))
    print(sys.argv[1], 'resid px', [round(v, 2) for v in d['resid']], 'xmax', round(d['xmax_mm']),
          [(c['L'], len(c['px']), round(c['mm'][-1][0]), round(c['mm'][-1][1])) for c in d['curves']])
