"""Independent verifier for the Farkas certificate written by certify_constant_gain_5_239.py.

It does not import the certifier.  From each row tag it rebuilds the row from the lemma
statements, checks the side conditions of the lemma with its own exact chord evaluator,
and then checks the Farkas combination exactly.

Usage: uv run --no-project python -I scripts/verify_constant_gain_5_239.py data/certificates/constant-gain-5-239-K160.json
"""
from fractions import Fraction as F
from math import isqrt
import bisect
import json
import sys

P12 = 10 ** 12


def rt_down(v):
    v = F(v)
    if v <= 0:
        return F(0)
    return F(isqrt(v.numerator * P12 * P12 // v.denominator), P12)


def rt_up(v):
    v = F(v)
    if v <= 0:
        return F(0)
    k = isqrt(v.numerator * P12 * P12 // v.denominator)
    while F(k, P12) ** 2 < v:
        k += 1
    return F(k, P12)


def cl(v):
    return F(-1) if v < -1 else (F(1) if v > 1 else v)


def chord_sq_box(a0, a1, r0, r1, s2, xw):
    """Squared same-side chord bound on [a0,a1] x [r0,r1], r0 > 0.
    Points q = p + r(cos phi, sin phi), |p| = a, phi in [0, pi], s2 <= |q|^2 <= 1:
    cos phi in [L, L + w],  L = (s2 - a^2 - r^2)/(2ar),  w = xw/(2ar).
    Arc angle f = arccos(cl L) - arccos(cl(L+w)) grows with w; on an L-interval its
    maximum is at an end or at L = -1 or L = 1 - w.  chord^2 = 2 r^2 (1 - cos f)."""
    Llo = (s2 - a1 * a1 - r1 * r1) / (2 * a0 * r0)
    Lhi = (s2 - a0 * a0 - r0 * r0) / (2 * a1 * r1)
    wmax = xw / (2 * a0 * r0)
    cands = {Llo, Lhi}
    if Llo <= -1 <= Lhi:
        cands.add(F(-1))
    if Llo <= 1 - wmax <= Lhi:
        cands.add(1 - wmax)
    worst_cos = F(1)
    for L in cands:
        u, v = cl(L), cl(L + wmax)
        cf = u * v + rt_down((1 - u * u) * (1 - v * v))
        worst_cos = min(worst_cos, cf)
    return 2 * r1 * r1 * (1 - worst_cos)


class Chords:
    def __init__(self, s, step=F(1, 2000), na=32):
        self.s = s
        s2 = s * s
        xw = 1 - s2
        self.edges = [step * k for k in range(int(2 / step) + 1)]
        agr = [s + (1 - s) * F(i, na) for i in range(na + 1)]
        self.tab = []
        for k in range(len(self.edges) - 1):
            r0, r1 = self.edges[k], self.edges[k + 1]
            if r0 == 0:
                self.tab.append(4 * r1 * r1)
                continue
            m = F(0)
            for i in range(na):
                m = max(m, chord_sq_box(agr[i], agr[i + 1], r0, r1, s2, xw))
            self.tab.append(m)

    def max_on(self, lo, hi):
        """max squared chord bound over table cells meeting (lo, hi] (coarse)"""
        i = max(bisect.bisect_right(self.edges, lo) - 1, 0)
        j = bisect.bisect_left(self.edges, hi)        # cell j-1 ends at edges[j] >= hi
        return max(self.tab[i:max(j, i + 1)])

    def _adapt(self, a0, a1, r0, r1, target, depth):
        if r1 <= 0:
            return True
        if r0 <= 0:
            return 4 * r1 * r1 <= target or (depth > 0 and
                   self._adapt(a0, a1, F(0), r1 / 2, target, depth - 1) and
                   self._adapt(a0, a1, r1 / 2, r1, target, depth - 1))
        if chord_sq_box(a0, a1, r0, r1, self.s * self.s, 1 - self.s * self.s) <= target:
            return True
        if depth == 0:
            return False
        am, rm = (a0 + a1) / 2, (r0 + r1) / 2
        return all(self._adapt(x0, x1, y0, y1, target, depth - 1)
                   for x0, x1 in ((a0, am), (am, a1)) for y0, y1 in ((r0, rm), (rm, r1)))

    def holds_on(self, lo, hi, target, na=32, depth=8):
        """True if every r in (lo, hi] has squared chord bound <= target (exact, adaptive)."""
        s = self.s
        agr = [s + (1 - s) * F(i, na) for i in range(na + 1)]
        i = max(bisect.bisect_right(self.edges, lo) - 1, 0)
        while i < len(self.tab) and self.edges[i] < hi:
            c0, c1 = self.edges[i], self.edges[i + 1]
            inside = (c0 >= lo and c1 <= hi)
            if not (inside and self.tab[i] <= target):
                r0, r1 = max(c0, lo), min(c1, hi)
                for k in range(na):
                    if not self._adapt(agr[k], agr[k + 1], r0, r1, target, depth):
                        return False
            i += 1
        return True


def main():
    cert = json.load(open(sys.argv[1]))
    delta = F(cert["delta"][0], cert["delta"][1])
    n0 = cert["n0"]
    G = [F(g) for g in cert["grid"]]
    assert G == sorted(G) and G[-1] == 2 and G[0] > 0
    eps = F(1, n0)
    xbar = 12 * delta + 4 * eps
    s = F(isqrt(int((1 - xbar) * 10 ** 12)), 10 ** 6)
    assert 0 < s and s * s <= 1 - xbar
    s2 = s * s
    xw = 1 - s2
    dbar = F(1, 3) + delta
    print(f"verifier: delta = {delta}, n0 = {n0}, xbar = {float(xbar):.6f}, s = {s}")
    ch = Chords(s)
    print(f"  chord table: {len(ch.tab)} cells, sup <= {float(rt_up(max(ch.tab))):.5f}")
    DD, EE = ('dd',), ('ee',)
    R = lambda g: ('R', g)
    A = lambda g: ('A', g)
    gset = set(G)

    def qpart(coef):   # coef * q,  q = 3dd - 1 + ee
        return {DD: 3 * coef, EE: coef}, -coef

    def merge(*ds):
        out = {}
        for d in ds:
            for k, v in d.items():
                out[k] = out.get(k, F(0)) + F(v)
        return out

    cvec = {}
    beta = F(0)
    nrows = 0
    for entry in cert["rows"]:
        tag = entry["tag"]
        y = F(entry["y"])
        assert y > 0
        fam = tag[0]
        p = [F(t) for t in tag[1:]]
        assert all(v in gset for v in p), ("tag parameter not a grid point", tag)
        if fam == 'Rmono':
            assert p[0] < p[1]
            d, rhs = {R(p[0]): 1, R(p[1]): -1}, F(0)
        elif fam == 'Amono':
            assert p[0] < p[1]
            d, rhs = {A(p[0]): 1, A(p[1]): -1}, F(0)
        elif fam == 'Xmono':
            assert p[0] < p[1]
            d, rhs = {R(p[0]): 3, A(p[0]): -1, R(p[1]): -3, A(p[1]): 1}, F(0)
        elif fam == 'Xpos':
            d, rhs = {A(p[0]): 1, R(p[0]): -3}, F(0)
        elif fam == 'Xleq':
            qd, k = qpart(-1)
            d, rhs = merge({R(p[0]): 3, A(p[0]): -1}, qd), -k
        elif fam == 'Rleqd':
            d, rhs = {R(p[0]): 1, DD: -1}, F(0)
        elif fam == 'R2':
            d, rhs = {R(F(2)): -1, DD: 1}, F(0)
        elif fam == 'A2lo':
            d, rhs = {A(F(2)): -1, EE: -1}, F(-1)
        elif fam == 'A2hi':
            d, rhs = {A(F(2)): 1, EE: 1}, F(1)
        elif fam in ('Pa', 'Pb', 'Pc'):
            rho, rp = p
            assert 0 < rho and rho / 2 < rp <= 2
            # Lemma C side condition: every r in (rho/2, r'] has chord <= rho
            assert ch.holds_on(rho / 2, rp, rho * rho), ("chord", tag)
            if rp == 2:
                base, rhs = {EE: -1, DD: -2}, F(-1)          # A(2) = 1 - ee, R(2) = dd
            else:
                base, rhs = {A(rp): 1, R(rp): -2}, F(0)
            d = merge({A(rho): -1}, base)
            if fam == 'Pb':
                d = merge(d, {A(rho / 2): 1, R(rho / 2): -1})
            if fam == 'Pc':
                d = merge(d, {A(rho / 2): -1, R(rho / 2): 2})
        elif fam == 'AA':
            lam, w, h = p
            assert 0 < lam < 2 and lam * lam >= 4 * xw and w > 0 and h <= 2
            r4 = rt_up(4 - lam * lam)
            assert w * w <= s2 + 1 - s * r4, ('w', tag)      # w^2 <= lo(lam)
            assert h * h >= 2 + r4, ('h', tag)               # h^2 >= hi(lam)
            assert ch.holds_on(w, F(2), lam * lam), ("chordA", tag)
            qd, k = qpart(-2)
            # separate dicts so that equal grid points (lam = h, lam = w, w = h) add up
            d = merge({EE: -2}, {A(lam): -2}, {A(h): F(-4, 3)}, {A(w): F(4, 3)}, qd)
            rhs = -2 - k
        elif fam == 'K':
            lo_, hi_ = p
            assert hi_ * hi_ <= 2 * s2 + 2 * s
            pts = [g for g in G if lo_ <= g <= hi_]
            d, rhs = {}, 2 * s
            for a_, b_ in zip(pts, pts[1:]):
                d[A(a_)] = d.get(A(a_), F(0)) + b_ * b_ - a_ * a_
        elif fam == 'Mlo':
            d = {}
            for a_, b_ in zip(G, G[1:]):
                w_ = (b_ * b_ - a_ * a_) / 2
                d = merge(d, {A(a_): -w_, R(a_): 6 * w_})
            qd, k = qpart(-4)
            d = merge(d, qd, {EE: 2})
            rhs = 1 - k
        elif fam == 'Mup':
            d = {}
            for a_, b_ in zip(G, G[1:]):
                w_ = (b_ * b_ - a_ * a_) / 2
                d = merge(d, {A(b_): w_, R(b_): -6 * w_})
            w0 = G[0] * G[0] / 2
            d = merge(d, {A(G[0]): w0, R(G[0]): -6 * w0})
            qd, k = qpart(4)
            d = merge(d, qd, {EE: -2})
            rhs = -1 - k
        else:
            raise ValueError(fam)
        for k_, v_ in d.items():
            cvec[k_] = cvec.get(k_, F(0)) + y * F(v_)
        beta += y * rhs
        nrows += 1
    # box: dd in [1/3 - 1/(3 n0), dbar], ee in [0, 1/n0], R in [0, dbar], A in [0, 1]
    mn = F(0)
    for k_, cj in cvec.items():
        if k_ == DD:
            lo, hi = F(1, 3) - eps / 3, dbar
        elif k_ == EE:
            lo, hi = F(0), eps
        elif k_[0] == 'R':
            lo, hi = F(0), dbar
        else:
            lo, hi = F(0), F(1)
        mn += cj * lo if cj > 0 else cj * hi
    gap = mn - beta
    print(f"  rows checked: {nrows}")
    print(f"  min_box (y^T A) v = {float(mn):.10f}, y^T b = {float(beta):.10f}, gap = {float(gap):.4e}")
    print("  VERIFIED: no profile satisfies the rows" if gap > 0 else "  FAILED")


if __name__ == "__main__":
    main()
