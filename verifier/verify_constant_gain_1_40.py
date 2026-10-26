"""Independent verifier for the case-split certificates written by certify_constant_gain_1_40.py.

It is scripts/verify_constant_gain_2_89.py extended by ONE row family (Lemma J, the
kappa-image length gap):
  J (lam, g2, gs):  A(g2) - A(lam) - A(gs) - 6 dd - 2 ee <= -2
side conditions (exact rationals): lam, g2, gs grid points; h0 = W/(2 lam) with h0 < s;
  g2 > lam;  g2^2 < 4 (s^2 - h0^2)  (so g2 < Lam0 = 2 sqrt(s^2 - h0^2));
  (gs + u)^2 >= 1 - h0^2 for a rational u <= sqrt(s^2 - h0^2), gs + u >= 0
  (so gs >= sig0 = sqrt(1 - h0^2) - sqrt(s^2 - h0^2)).
Everything else is the banked verifier verbatim (S2 families + Lemma AB + case split).
It does not import any certifier; it imports only the banked chord evaluator.

Usage: uv run --no-project python -I scripts/verify_constant_gain_1_40.py data/certificates/constant-gain-1-40-K160-c60.json.gz
"""
from fractions import Fraction as F
from math import isqrt
import gzip
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_constant_gain_5_239 import rt_up, Chords      # banked chord evaluator (no certifier code)


def sqrt_floor_rat(q, prec=10 ** 12):
    """rational u with 0 <= u <= sqrt(q) (q >= 0 rational)."""
    q = F(q)
    if q <= 0:
        return F(0)
    u = F(isqrt(q.numerator * prec * prec // q.denominator), prec)
    assert u * u <= q
    return u

DD, EE = ('dd',), ('ee',)


def R(g):
    return ('R', g)


def A(g):
    return ('A', g)


def add_into(d, other, c=F(1)):
    for k, v in other.items():
        d[k] = d.get(k, F(0)) + c * v


def main():
    path = sys.argv[1]
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        cert = json.load(fh)
    delta = F(cert["delta"][0], cert["delta"][1])
    n0 = cert["n0"]
    G = [F(g) for g in cert["grid"]]
    assert G == sorted(G) and G[-1] == 2 and G[0] > 0 and len(set(G)) == len(G)
    gset = set(G)
    eps = F(1, n0)
    xbar = 12 * delta + 4 * eps
    s = F(isqrt(int((1 - xbar) * 10 ** 12)), 10 ** 6)
    assert 0 < s and s * s <= 1 - xbar
    s2 = s * s
    xw = 1 - s2                     # = W
    dbar = F(1, 3) + delta
    print(f"verifier: delta = {delta}, n0 = {n0}, xbar = {float(xbar):.6f}, s = {s}, W = {float(xw):.6f}")
    ch = Chords(s)
    print(f"  chord table: {len(ch.tab)} cells, sup <= {float(rt_up(max(ch.tab))):.5f}")
    lams = [F(l) for l in cert["lams"]]
    edges = {F(l): [F(e) for e in es] for l, es in cert["edges"].items()}
    assert set(edges) == set(lams), "every AB threshold must pass the side-condition check"
    # ---- side conditions of Lemma AB, once per threshold
    for lam in lams:
        assert lam in gset and 0 < lam < 2, lam
        assert ch.holds_on(F(0), F(2), lam * lam), ("t_all > lam", lam)   # every 3-point circle has a pair <= lam
        E = edges[lam]
        assert all(e in gset for e in E) and E == sorted(E) and len(set(E)) == len(E)
        assert E[0] <= lam / 2 and E[-1] == 2
    print(f"  Lemma AB side conditions: {len(lams)} thresholds OK (t_all <= lam, bands cover (e0, 2], e0 <= lam/2)")

    def qpart(coef):
        return {DD: 3 * coef, EE: coef}, -coef

    up_cache = {}

    def int_rhoX(upper):
        if upper in up_cache:
            return up_cache[upper]
        d = {}
        for a_, b_ in zip(G, G[1:]):
            w = (b_ * b_ - a_ * a_) / 2
            g = b_ if upper else a_
            d[R(g)] = d.get(R(g), F(0)) + 3 * w
            d[A(g)] = d.get(A(g), F(0)) - w
        if upper:
            w0 = G[0] * G[0] / 2
            d[R(G[0])] = d.get(R(G[0]), F(0)) + 3 * w0
            d[A(G[0])] = d.get(A(G[0]), F(0)) - w0
        up_cache[upper] = d
        return d

    def build(tag, case):
        fam = tag[0]
        if fam in ('AB1', 'AB2', 'AB3', 'AB4', 'AB5u', 'AB5v', 'AB6u', 'AB6v', 'CASElo', 'CASEhi'):
            if fam == 'CASElo':
                Ylo = F(tag[1]); assert Ylo == case[0]
                # Ylo <= W - 4q + 2 Up
                qd, k = qpart(F(4)); d = dict(qd); add_into(d, int_rhoX(True), F(-2))
                return d, xw - Ylo - k
            if fam == 'CASEhi':
                Yhi = F(tag[1]); assert Yhi == case[1]
                qd, k = qpart(F(-4)); d = dict(qd); add_into(d, int_rhoX(False), F(2))
                return d, Yhi - xw - k
            lam = F(tag[1]); assert lam in edges
            E = edges[lam]; m = len(E) - 1
            if fam == 'AB1':
                qd, k = qpart(F(-2))
                d = {A(lam): F(-2), EE: F(-2)}; add_into(d, qd)
                for j in range(m):
                    d[('W', lam, j)] = F(-1)
                return d, -2 - k
            if fam in ('AB5u', 'AB5v'):
                key = 'U' if fam == 'AB5u' else 'V'
                coef = F(4) if fam == 'AB5u' else F(4, 3)
                d = {(key, lam, j): F(1) for j in range(m)}
                qd, k = qpart(4 * coef); add_into(d, qd); add_into(d, int_rhoX(True), -2 * coef)
                return d, coef * xw - k
            if fam == 'AB6u':
                Yhi = F(tag[2]); assert Yhi == case[1]
                d = {('U', lam, j): F(1) for j in range(m)}
                d[DD] = -12 * Yhi; d[R(lam)] = 12 * Yhi
                return d, F(0)
            if fam == 'AB6v':
                h = F(tag[2]); Yhi = F(tag[3]); assert Yhi == case[1]
                assert h in gset and h <= lam / 2
                d = {('V', lam, j): F(1) for j in range(m)}
                d[DD] = -4 * Yhi; d[R(h)] = 4 * Yhi
                return d, F(0)
            j = int(tag[2]); assert 0 <= j < m
            ea, eb = E[j], E[j + 1]
            if fam == 'AB2':
                d = {('W', lam, j): F(1), A(eb): F(-4, 3)}
                d[A(ea)] = d.get(A(ea), F(0)) + F(4, 3)
                return d, F(0)
            t = F(tag[3]); assert t > 0
            c0 = 4 * s2 - lam * lam
            if fam == 'AB3':
                g = 2 * ea * ea - ((4 + t) * s2 + c0 / t)
                return {('W', lam, j): g, ('U', lam, j): -(1 + 2 / t), ('V', lam, j): -(2 + t)}, F(0)
            if fam == 'AB4':
                g = ((4 - t) * s2 - c0 / t) - 2 * eb * eb
                cu = max(2 / t - 1, F(0)); cv = max(t - 2, F(0))
                return {('W', lam, j): g, ('U', lam, j): -cu, ('V', lam, j): -cv}, F(0)
        # ---- banked families (verbatim logic of scripts/verify_constant_gain_5_239.py)
        p = [F(x) for x in tag[1:]]
        assert all(v in gset for v in p), ("tag parameter not a grid point", tag)
        if fam == 'Rmono':
            assert p[0] < p[1]; return {R(p[0]): F(1), R(p[1]): F(-1)}, F(0)
        if fam == 'Amono':
            assert p[0] < p[1]; return {A(p[0]): F(1), A(p[1]): F(-1)}, F(0)
        if fam == 'Xmono':
            assert p[0] < p[1]; return {R(p[0]): F(3), A(p[0]): F(-1), R(p[1]): F(-3), A(p[1]): F(1)}, F(0)
        if fam == 'Xpos':
            return {A(p[0]): F(1), R(p[0]): F(-3)}, F(0)
        if fam == 'Xleq':
            qd, k = qpart(F(-1)); d = {R(p[0]): F(3), A(p[0]): F(-1)}; add_into(d, qd); return d, -k
        if fam == 'Rleqd':
            return {R(p[0]): F(1), DD: F(-1)}, F(0)
        if fam == 'R2':
            return {R(F(2)): F(-1), DD: F(1)}, F(0)
        if fam == 'A2lo':
            return {A(F(2)): F(-1), EE: F(-1)}, F(-1)
        if fam == 'A2hi':
            return {A(F(2)): F(1), EE: F(1)}, F(1)
        if fam in ('Pa', 'Pb', 'Pc'):
            rho, rp = p
            assert 0 < rho and rho / 2 < rp <= 2
            assert ch.holds_on(rho / 2, rp, rho * rho), ("chord", tag)
            if rp == 2:
                d, rhs = {EE: F(-1), DD: F(-2)}, F(-1)
            else:
                d, rhs = {A(rp): F(1), R(rp): F(-2)}, F(0)
            d[A(rho)] = d.get(A(rho), F(0)) - 1
            if fam == 'Pb':
                add_into(d, {A(rho / 2): F(1), R(rho / 2): F(-1)})
            if fam == 'Pc':
                add_into(d, {A(rho / 2): F(-1), R(rho / 2): F(2)})
            return d, rhs
        if fam == 'AA':
            lam, w, h = p
            assert 0 < lam < 2 and lam * lam >= 4 * xw and w > 0 and h <= 2
            r4 = rt_up(4 - lam * lam)
            assert w * w <= s2 + 1 - s * r4 and h * h >= 2 + r4
            assert ch.holds_on(w, F(2), lam * lam)
            qd, k = qpart(F(-2))
            d = {EE: F(-2), A(lam): F(-2)}
            d[A(h)] = d.get(A(h), F(0)) - F(4, 3); d[A(w)] = d.get(A(w), F(0)) + F(4, 3)
            add_into(d, qd)
            return d, -2 - k
        if fam == 'K':
            lo_, hi_ = p
            assert hi_ * hi_ <= 2 * s2 + 2 * s
            pts = [g for g in G if lo_ <= g <= hi_]
            d = {}
            for a_, b_ in zip(pts, pts[1:]):
                d[A(a_)] = d.get(A(a_), F(0)) + b_ * b_ - a_ * a_
            return d, 2 * s
        if fam == 'Mlo':
            d = {}
            for a_, b_ in zip(G, G[1:]):
                w_ = (b_ * b_ - a_ * a_) / 2
                d[A(a_)] = d.get(A(a_), F(0)) - w_; d[R(a_)] = d.get(R(a_), F(0)) + 6 * w_
            qd, k = qpart(F(-4)); add_into(d, qd); d[EE] = d.get(EE, F(0)) + 2
            return d, 1 - k
        if fam == 'Mup':
            d = {}
            for a_, b_ in zip(G, G[1:]):
                w_ = (b_ * b_ - a_ * a_) / 2
                d[A(b_)] = d.get(A(b_), F(0)) + w_; d[R(b_)] = d.get(R(b_), F(0)) - 6 * w_
            w0 = G[0] * G[0] / 2
            d[A(G[0])] = d.get(A(G[0]), F(0)) + w0; d[R(G[0])] = d.get(R(G[0]), F(0)) - 6 * w0
            qd, k = qpart(F(4)); add_into(d, qd); d[EE] = d.get(EE, F(0)) - 2
            return d, -1 - k
        if fam == 'J':
            lam, g2, gs = p
            h0 = xw / (2 * lam)
            assert 0 <= h0 < s, ("J: h0 >= s", tag)
            assert g2 > lam, ("J: g2 <= lam", tag)
            assert g2 * g2 < 4 * (s2 - h0 * h0), ("J: g2 >= Lam0", tag)
            u = sqrt_floor_rat(s2 - h0 * h0)
            assert gs + u >= 0 and (gs + u) * (gs + u) >= 1 - h0 * h0, ("J: gs < sig0", tag)
            d = {}
            for k_, v_ in ((A(g2), F(1)), (A(lam), F(-1)), (A(gs), F(-1)), (DD, F(-6)), (EE, F(-2))):
                d[k_] = d.get(k_, F(0)) + v_
            return d, F(-2)
        raise ValueError(fam)

    # ---- cases cover [0, W]
    cases = [(F(c["Ylo"]), F(c["Yhi"])) for c in cert["cases"]]
    assert cases[0][0] == 0 and cases[-1][1] == xw
    for (a0, b0), (a1, b1) in zip(cases, cases[1:]):
        assert a0 < b0 and b0 == a1
    print(f"  {len(cases)} cases cover [0, W]")
    worst = None
    for c, case in zip(cert["cases"], cases):
        cvec = {}
        beta = F(0)
        for entry in c["rows"]:
            y = F(entry["y"]); assert y > 0
            d, rhs = build(entry["tag"], case)
            for k_, v_ in d.items():
                cvec[k_] = cvec.get(k_, F(0)) + y * v_
            beta += y * rhs
        mn = F(0)
        for k_, cj in cvec.items():
            if k_ == DD:
                lo, hi = F(1, 3) - eps / 3, dbar
            elif k_ == EE:
                lo, hi = F(0), eps
            elif k_[0] == 'R':
                lo, hi = F(0), dbar
            elif k_[0] == 'A':
                lo, hi = F(0), F(1)
            elif k_[0] == 'W':
                lo, hi = F(0), F(2)
            elif k_[0] in ('U', 'V'):
                lo, hi = F(0), F(4)
            else:
                raise ValueError(k_)
            mn += cj * lo if cj > 0 else cj * hi
        gap = mn - beta
        print(f"  case [{float(case[0]):.4f},{float(case[1]):.4f}]: rows {len(c['rows'])}, gap = {float(gap):.4e} "
              f"{'ok' if gap > 0 else 'FAIL'}", flush=True)
        if gap <= 0:
            print("  FAILED"); return
        worst = gap if worst is None or gap < worst else worst
    print(f"  smallest gap = {float(worst):.4e}")
    print("  VERIFIED: in every case no profile satisfies the rows")


if __name__ == "__main__":
    main()
