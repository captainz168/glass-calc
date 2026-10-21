"""Nonlinear (von Karman) rectangular plate, Ritz (Legendre) + GFPM risk.
Non-dimensional: t = 1, E = 1. Plate x in [-a/2, a/2], y in [-b/2, b/2].
Supported edges simple (w = 0), all edges free in-plane.
edges: subset of 'LRBT' (L: x=-a/2, R: x=+a/2, B: y=-b/2, T: y=+b/2).
Symmetry is exploited automatically: x-symmetric if L,R same; y-symmetric if B,T same.
"""
import numpy as np
from numpy.polynomial import legendre as Lg

def leg(n, xi):
    V = np.zeros((n + 1, len(xi))); D1 = np.zeros_like(V); D2 = np.zeros_like(V)
    for i in range(n + 1):
        c = np.zeros(i + 1); c[i] = 1
        V[i] = Lg.legval(xi, c); D1[i] = Lg.legval(xi, Lg.legder(c, 1)); D2[i] = Lg.legval(xi, Lg.legder(c, 2))
    return V, D1, D2

def basis1d(n, xi, lo, hi, parity):
    """f_i = g * P_i, g=(1+xi)^lo (1-xi)^hi ; parity: None/'even'/'odd' selects P_i."""
    V, D1, D2 = leg(n, xi)
    lo = int(bool(lo)); hi = int(bool(hi))
    g = (1 + xi) ** lo * (1 - xi) ** hi
    if lo and hi: g1, g2 = -2 * xi, -2 * np.ones_like(xi)
    elif lo: g1, g2 = np.ones_like(xi), np.zeros_like(xi)
    elif hi: g1, g2 = -np.ones_like(xi), np.zeros_like(xi)
    else: g1, g2 = np.zeros_like(xi), np.zeros_like(xi)
    F = g * V; F1 = g1 * V + g * D1; F2 = g2 * V + 2 * g1 * D1 + g * D2
    idx = list(range(n + 1))
    if parity == 'even': idx = [i for i in idx if i % 2 == 0]
    if parity == 'odd': idx = [i for i in idx if i % 2 == 1]
    return F[idx], F1[idx], F2[idx], idx

class Plate:
    def __init__(self, a, b, edges, nw=20, nuv=18, ng=40, nu=0.22):
        self.a, self.b, self.nu, self.edges = a, b, nu, edges
        symx = (('L' in edges) == ('R' in edges)); symy = (('B' in edges) == ('T' in edges))
        self.symx, self.symy = symx, symy
        xg, wg = Lg.leggauss(2 * ng)
        # quadrature: if symmetric, integrate half [0,1] and multiply by 2
        if symx: mx = xg >= 0; xq, wx = xg[mx], wg[mx] * 2
        else: xq, wx = xg, wg
        if symy: my = xg >= 0; yq, wy = xg[my], wg[my] * 2
        else: yq, wy = xg, wg
        self.xq, self.yq = xq, yq
        sx, sy = 2 / a, 2 / b
        self.dA = (np.outer(wx, wy) * (a / 2) * (b / 2)).ravel()
        px = lambda even: (('even' if even else 'odd') if symx else None)
        py = lambda even: (('even' if even else 'odd') if symy else None)
        T = lambda A, B: np.einsum('ip,jq->pqij', A, B).reshape(len(xq) * len(yq), -1)
        fx, fx1, fx2, _ = basis1d(nw, xq, 'L' in edges, 'R' in edges, px(True))
        fy, fy1, fy2, _ = basis1d(nw, yq, 'B' in edges, 'T' in edges, py(True))
        self.Bw = T(fx, fy); self.Bwx = T(fx1, fy) * sx; self.Bwy = T(fx, fy1) * sy
        self.Bwxx = T(fx2, fy) * sx * sx; self.Bwyy = T(fx, fy2) * sy * sy; self.Bwxy = T(fx1, fy1) * sx * sy
        # u: odd in x (if symx), even in y ; v: even in x, odd in y
        ux, ux1, _, ix = basis1d(nuv, xq, 0, 0, px(False)); uy, uy1, _, iy = basis1d(nuv, yq, 0, 0, py(True))
        Bu, Bux, Buy = T(ux, uy), T(ux1, uy) * sx, T(ux, uy1) * sy
        keep = np.ones(Bu.shape[1], bool)
        if not symx: keep[0] = False  # remove rigid translation (P0*P0)
        self.Bu, self.Bux, self.Buy = Bu[:, keep], Bux[:, keep], Buy[:, keep]
        vx, vx1, _, _ = basis1d(nuv, xq, 0, 0, px(True)); vy, vy1, _, _ = basis1d(nuv, yq, 0, 0, py(False))
        Bv, Bvx, Bvy = T(vx, vy), T(vx1, vy) * sx, T(vx, vy1) * sy
        keep = np.ones(Bv.shape[1], bool)
        if not symy: keep[0] = False
        self.Bv, self.Bvx, self.Bvy = Bv[:, keep], Bvx[:, keep], Bvy[:, keep]
        self.rot = not (symx or symy)
        self.nw, self.nu_, self.nv = self.Bw.shape[1], self.Bu.shape[1], self.Bv.shape[1]
        self.n = self.nw + self.nu_ + self.nv
        D = 1 / (12 * (1 - nu ** 2)); self.D = D; self.Cm = 1 / (1 - nu ** 2)
        dA = self.dA[:, None]
        self.Kb = D * ((self.Bwxx * dA).T @ self.Bwxx + (self.Bwyy * dA).T @ self.Bwyy
                       + nu * ((self.Bwxx * dA).T @ self.Bwyy + (self.Bwyy * dA).T @ self.Bwxx)
                       + 2 * (1 - nu) * (self.Bwxy * dA).T @ self.Bwxy)
        self.Fq = (self.Bw * dA).sum(0)
        self.c = np.zeros(self.n)

    def split(self, c):
        return c[:self.nw], c[self.nw:self.nw + self.nu_], c[self.nw + self.nu_:]

    def strains(self, c):
        W, U, V = self.split(c)
        wx = self.Bwx @ W; wy = self.Bwy @ W
        ex = self.Bux @ U + 0.5 * wx ** 2; ey = self.Bvy @ V + 0.5 * wy ** 2
        g = self.Buy @ U + self.Bvx @ V + wx * wy
        return wx, wy, ex, ey, g

    def gradhess(self, c, q, hess=True):
        W = c[:self.nw]
        wx, wy, ex, ey, gm = self.strains(c)
        Cm, nu, dA = self.Cm, self.nu, self.dA
        Nx = Cm * (ex + nu * ey); Ny = Cm * (ey + nu * ex); Nxy = Cm * (1 - nu) / 2 * gm
        Zu = np.zeros((len(dA), self.nu_)); Zv = np.zeros((len(dA), self.nv))
        Jex = np.hstack([wx[:, None] * self.Bwx, self.Bux, Zv])
        Jey = np.hstack([wy[:, None] * self.Bwy, Zu, self.Bvy])
        Jg = np.hstack([wy[:, None] * self.Bwx + wx[:, None] * self.Bwy, self.Buy, self.Bvx])
        g = Jex.T @ (Nx * dA) + Jey.T @ (Ny * dA) + Jg.T @ (Nxy * dA)
        g[:self.nw] += self.Kb @ W - q * self.Fq
        if not hess: return g, None
        Gm = (Cm * dA)[:, None]
        H = (Jex * Gm).T @ (Jex + nu * Jey) + (Jey * Gm).T @ (Jey + nu * Jex) + (Jg * (Gm * (1 - nu) / 2)).T @ Jg
        nw = self.nw
        H[:nw, :nw] += self.Kb + (self.Bwx * (Nx * dA)[:, None]).T @ self.Bwx + (self.Bwy * (Ny * dA)[:, None]).T @ self.Bwy \
            + (self.Bwx * (Nxy * dA)[:, None]).T @ self.Bwy + (self.Bwy * (Nxy * dA)[:, None]).T @ self.Bwx
        if self.rot:
            H[nw:, nw:] += np.eye(self.n - nw) * 1e-10 * self.a * self.b
        return g, H

    def _newton(self, c, q, tol, maxit):
        ref = abs(q) * np.linalg.norm(self.Fq) + 1e-300
        for it in range(maxit):
            g, H = self.gradhess(c, q)
            r = np.linalg.norm(g) / ref
            if r < tol: return c, True, it
            if not np.isfinite(r) or r > 1e8: return c, False, it
            c = c + np.linalg.solve(H, -g)
        return c, False, maxit

    def solve(self, q, tol=1e-9, maxit=25):
        """Newton with automatic load sub-stepping from current state (self.q)."""
        q0 = getattr(self, 'q', 0.0); c0 = self.c.copy()
        steps = [q]
        while steps:
            qt = steps[0]
            c, ok, it = self._newton(c0.copy(), qt, tol, maxit)
            if ok:
                c0, q0 = c, qt; steps.pop(0)
            else:
                steps.insert(0, (q0 + qt) / 2)
                if len(steps) > 40: raise RuntimeError('no convergence')
        self.c, self.q = c0, q0
        return it

    def face_stresses(self):
        W = self.c[:self.nw]; nu = self.nu
        wx, wy, ex, ey, gm = self.strains(self.c)
        k = 1 / (1 - nu ** 2)
        smx = k * (ex + nu * ey); smy = k * (ey + nu * ex); smxy = gm / (2 * (1 + nu))
        kxx = self.Bwxx @ W; kyy = self.Bwyy @ W; kxy = self.Bwxy @ W
        z = 0.5
        sbx = -k * z * (kxx + nu * kyy); sby = -k * z * (kyy + nu * kxx); sbxy = -z * kxy / (1 + nu)
        out = []
        for s in (+1, -1):
            a_ = smx + s * sbx; b_ = smy + s * sby; c_ = smxy + s * sbxy
            cen = (a_ + b_) / 2; rad = np.sqrt(((a_ - b_) / 2) ** 2 + c_ ** 2)
            out.append((cen + rad, cen - rad))
        return out

    def wmax(self):
        # center / max deflection on a fine grid of quadrature points
        return np.max(np.abs(self.Bw @ self.c[:self.nw]))

def cfac(r):
    r = np.clip(r, -1, 1)
    return -0.005 * r ** 6 + 0.022 * r ** 5 + 0.055 * r ** 4 + 0.039 * r ** 3 + 0.031 * r ** 2 + 0.06 * r + 0.8

def risk_integral(p, m=7):
    """returns S = int (c*s1)^m dA over each face (tuple)"""
    res = []
    for s1, s2 in p.face_stresses():
        pos = s1 > 0
        r = np.where(pos, s2 / np.where(pos, s1, 1.0), 0.0)
        res.append(np.sum(np.where(pos, (cfac(r) * s1) ** m, 0.0) * p.dA))
    return res
