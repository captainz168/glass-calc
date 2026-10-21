import json, numpy as np
E = 71.7e9; K = 2.86e-53; M = 7; DUR = (3 / 60) ** (7 / 16)
TARGET = -np.log(1 - 0.008)

class Table:
    def __init__(self, fn, faces='both'):
        d = json.load(open(fn))
        self.lq = np.log(np.array(d['qhat'])); self.AR = np.array(d['AR'])
        J1 = np.array(d['J']); J2 = np.array(d['J2'])
        self.J = np.log(np.exp(J1) + np.exp(J2)) if faces == 'both' else J1
        self.w = np.array(d['w'])

    def _row(self, arr, AR):
        x = np.log(self.AR); xa = np.log(np.clip(AR, self.AR[0], self.AR[-1]))
        i = np.clip(np.searchsorted(x, xa) - 1, 0, len(x) - 2)
        f = (xa - x[i]) / (x[i + 1] - x[i])
        return arr[i] * (1 - f) + arr[i + 1] * f

    def nfl(self, a_mm, b_mm, t_mm, AR=None):
        """a,b: plate x/y lengths (mm). AR as used in table (default b/a)."""
        ab = a_mm * b_mm * 1e-6; t = t_mm * 1e-3
        AR = b_mm / a_mm if AR is None else AR
        J = self._row(self.J, AR)
        need = np.log(TARGET / (K * DUR)) - (1 - M) * np.log(ab) - M * np.log(E * t * t)
        lqh = np.interp(need, J, self.lq)
        q = np.exp(lqh) * E * t ** 4 / ab ** 2
        return q / 1e3  # kPa

    def defl(self, a_mm, b_mm, t_mm, q_kpa, AR=None):
        ab = a_mm * b_mm * 1e-6; t = t_mm * 1e-3
        AR = b_mm / a_mm if AR is None else AR
        qh = q_kpa * 1e3 * ab ** 2 / (E * t ** 4)
        return np.interp(np.log(qh), self.lq, self._row(self.w, AR)) * t_mm
