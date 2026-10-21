"""Two-ply laminated glass plate (glass/PVB/glass), von Karman, Ritz.
Common w; each ply has own in-plane displacements; interlayer carries shear only.
Units: SI (m, Pa). 4-sided simple support, all edges free in-plane, quarter symmetry.
"""
import numpy as np
from numpy.polynomial import legendre as Lg
from plate import basis1d, cfac

class LamPlate:
    def __init__(self, a, b, h1, h2, hi, G, E=71.7e9, nu=0.22, edges='LRBT', nw=20, nuv=18, ng=32):
        self.a, self.b, self.h = a, b, (h1, h2)
        self.hi, self.G, self.E, self.nu = hi, G, E, nu
        symx = (('L' in edges) == ('R' in edges)); symy = (('B' in edges) == ('T' in edges))
        assert symx and symy
        xg, wg = Lg.leggauss(2 * ng); mk = xg >= 0
        xq, wq = xg[mk], wg[mk] * 2
        self.dA = (np.outer(wq, wq) * (a / 2) * (b / 2)).ravel()
        sx, sy = 2 / a, 2 / b
        T = lambda A, B: np.einsum('ip,jq->pqij', A, B).reshape(len(xq) * len(xq), -1)
        fx, fx1, fx2, _ = basis1d(nw, xq, 'L' in edges, 'R' in edges, 'even')
        fy, fy1, fy2, _ = basis1d(nw, xq, 'B' in edges, 'T' in edges, 'even')
        self.Bw = T(fx, fy); self.Bwx = T(fx1, fy) * sx; self.Bwy = T(fx, fy1) * sy
        self.Bwxx = T(fx2, fy) * sx * sx; self.Bwyy = T(fx, fy2) * sy * sy; self.Bwxy = T(fx1, fy1) * sx * sy
        ux, ux1, _, _ = basis1d(nuv, xq, 0, 0, 'odd'); uy, uy1, _, _ = basis1d(nuv, xq, 0, 0, 'even')
        self.Bu, self.Bux, self.Buy = T(ux, uy), T(ux1, uy) * sx, T(ux, uy1) * sy
        vx, vx1, _, _ = basis1d(nuv, xq, 0, 0, 'even'); vy, vy1, _, _ = basis1d(nuv, xq, 0, 0, 'odd')
        self.Bv, self.Bvx, self.Bvy = T(vx, vy), T(vx1, vy) * sx, T(vx, vy1) * sy
        self.nw = self.Bw.shape[1]; self.nm = self.Bu.shape[1]
        # unknown layout: W, U1, V1, U2, V2
        self.sl = {}
        o = 0
        for k, n in [('W', self.nw), ('U1', self.nm), ('V1', self.nm), ('U2', self.nm), ('V2', self.nm)]:
            self.sl[k] = slice(o, o + n); o += n
        self.n = o
        dA = self.dA[:, None]
        Kb0 = ((self.Bwxx * dA).T @ self.Bwxx + (self.Bwyy * dA).T @ self.Bwyy
               + nu * ((self.Bwxx * dA).T @ self.Bwyy + (self.Bwyy * dA).T @ self.Bwxx)
               + 2 * (1 - nu) * (self.Bwxy * dA).T @ self.Bwxy)
        Dsum = sum(E * hk ** 3 / (12 * (1 - nu ** 2)) for hk in (h1, h2))
        self.Kb = Dsum * Kb0
        self.Fq = (self.Bw * dA).sum(0)
        self.d = hi + (h1 + h2) / 2
        self.c = np.zeros(self.n); self.q = 0.0

    def _ply(self, c, k):
        W = c[self.sl['W']]; U = c[self.sl['U%d' % k]]; V = c[self.sl['V%d' % k]]
        wx = self.Bwx @ W; wy = self.Bwy @ W
        ex = self.Bux @ U + 0.5 * wx ** 2; ey = self.Bvy @ V + 0.5 * wy ** 2
        g = self.Buy @ U + self.Bvx @ V + wx * wy
        return wx, wy, ex, ey, g

    def gradhess(self, c, q, hess=True):
        E, nu, dA, n = self.E, self.nu, self.dA, self.n
        g = np.zeros(n); H = np.zeros((n, n)) if hess else None
        sW = self.sl['W']; nw = self.nw
        for k in (1, 2):
            hk = self.h[k - 1]; Cm = E * hk / (1 - nu ** 2)
            wx, wy, ex, ey, gm = self._ply(c, k)
            Nx = Cm * (ex + nu * ey); Ny = Cm * (ey + nu * ex); Nxy = Cm * (1 - nu) / 2 * gm
            Jex = np.zeros((len(dA), n)); Jey = np.zeros_like(Jex); Jg = np.zeros_like(Jex)
            Jex[:, sW] = wx[:, None] * self.Bwx; Jex[:, self.sl['U%d' % k]] = self.Bux
            Jey[:, sW] = wy[:, None] * self.Bwy; Jey[:, self.sl['V%d' % k]] = self.Bvy
            Jg[:, sW] = wy[:, None] * self.Bwx + wx[:, None] * self.Bwy
            Jg[:, self.sl['U%d' % k]] = self.Buy; Jg[:, self.sl['V%d' % k]] = self.Bvx
            g += Jex.T @ (Nx * dA) + Jey.T @ (Ny * dA) + Jg.T @ (Nxy * dA)
            if hess:
                Gm = (Cm * dA)[:, None]
                H += (Jex * Gm).T @ (Jex + nu * Jey) + (Jey * Gm).T @ (Jey + nu * Jex) + (Jg * (Gm * (1 - nu) / 2)).T @ Jg
                H[sW, sW] += (self.Bwx * (Nx * dA)[:, None]).T @ self.Bwx + (self.Bwy * (Ny * dA)[:, None]).T @ self.Bwy \
                    + (self.Bwx * (Nxy * dA)[:, None]).T @ self.Bwy + (self.Bwy * (Nxy * dA)[:, None]).T @ self.Bwx
        # interlayer shear: gx = (u1 - u2 + d wx)/hi, gy = (v1 - v2 + d wy)/hi ; energy 0.5*G*hi*(gx^2+gy^2)
        W = c[sW]
        Gk = self.G / self.hi
        Jx = np.zeros((len(dA), n)); Jy = np.zeros_like(Jx)
        Jx[:, sW] = self.d * self.Bwx; Jx[:, self.sl['U1']] = self.Bu; Jx[:, self.sl['U2']] = -self.Bu
        Jy[:, sW] = self.d * self.Bwy; Jy[:, self.sl['V1']] = self.Bv; Jy[:, self.sl['V2']] = -self.Bv
        gx = Jx @ c; gy = Jy @ c
        g += Gk * (Jx.T @ (gx * dA) + Jy.T @ (gy * dA))
        g[sW] += self.Kb @ W - q * self.Fq
        if hess:
            H += Gk * ((Jx * dA[:, None]).T @ Jx + (Jy * dA[:, None]).T @ Jy)
            H[sW, sW] += self.Kb
        return g, H

    def _newton(self, c, q, tol=1e-9, maxit=25):
        ref = abs(q) * np.linalg.norm(self.Fq) + 1e-300
        for it in range(maxit):
            g, H = self.gradhess(c, q)
            r = np.linalg.norm(g) / ref
            if r < tol: return c, True
            if not np.isfinite(r) or r > 1e8: return c, False
            c = c + np.linalg.solve(H, -g)
        return c, False

    def solve(self, q):
        q0, c0 = self.q, self.c.copy(); steps = [q]
        while steps:
            c, ok = self._newton(c0.copy(), steps[0])
            if ok: c0, q0 = c, steps.pop(0)
            else:
                steps.insert(0, (q0 + steps[0]) / 2)
                if len(steps) > 40: raise RuntimeError('no conv')
        self.c, self.q = c0, q0

    def surfaces(self):
        """principal stresses (s1,s2) for the 4 glass surfaces"""
        E, nu = self.E, self.nu; W = self.c[self.sl['W']]
        kxx = self.Bwxx @ W; kyy = self.Bwyy @ W; kxy = self.Bwxy @ W
        k_ = E / (1 - nu ** 2); out = []
        for k in (1, 2):
            hk = self.h[k - 1]
            wx, wy, ex, ey, gm = self._ply(self.c, k)
            smx = k_ * (ex + nu * ey); smy = k_ * (ey + nu * ex); smxy = E / (2 * (1 + nu)) * gm
            for s in (+1, -1):
                z = s * hk / 2
                ax = smx - k_ * z * (kxx + nu * kyy); ay = smy - k_ * z * (kyy + nu * kxx); axy = smxy - E / (1 + nu) * z * kxy
                cen = (ax + ay) / 2; rad = np.sqrt(((ax - ay) / 2) ** 2 + axy ** 2)
                out.append((cen + rad, cen - rad))
        return out

    def risk(self, m=7):
        res = []
        for s1, s2 in self.surfaces():
            pos = s1 > 0
            r = np.where(pos, s2 / np.where(pos, s1, 1.0), 0.0)
            res.append(np.sum(np.where(pos, (cfac(r) * s1) ** m, 0.0) * self.dA))
        return res

    def wmax(self):
        return np.max(np.abs(self.Bw @ self.c[self.sl['W']]))
