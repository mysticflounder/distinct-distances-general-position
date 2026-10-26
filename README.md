# Distinct distances in general position: more than (1/3 + 1/40)n

**Theorem.** Let P be a set of n ≥ 10⁴ points in the plane, with no three
points collinear and no four points concyclic. Then P determines more than
(1/3 + 1/40)n distinct distances.

Every circle centred at a point of such a set holds at most three other
points, so the set determines at least (n − 1)/3 distinct distances. The
theorem improves this coefficient 1/3 by 1/40.

The proof normalizes the set into a thin annulus, derives linear inequalities
for the distribution of the pairwise distances, and closes with an exact
rational Farkas certificate over 60 cases. An independent exact verifier
checks the certificate. The Lean 4 kernel checks the complete proof, with the
certificate.

**Status:** initial draft of the paper (October 2026). The page ranges in the
bibliography are not yet checked against the primary sources.

## Contents

| Path | Contents |
|---|---|
| `paper/distinct-distances-general-position.tex` | the paper (LaTeX) |
| `certificate/constant-gain-1-40-K160-c60.json.gz` | the certificate: 60 cases, rational multipliers |
| `verifier/verify_constant_gain_1_40.py` | exact verifier (Python standard library only) |
| `verifier/verify_constant_gain_5_239.py` | chord evaluator that the verifier imports |
| `lean/` | Lean 4 proof of the theorem |

## Check the certificate

```sh
python3 -I verifier/verify_constant_gain_1_40.py certificate/constant-gain-1-40-K160-c60.json.gz
```

The verifier rebuilds every row from its tag, checks every side condition,
and evaluates each case in exact rational arithmetic. It prints `VERIFIED` on
success. With Python 3.13 it uses about 92 s of CPU time on one core and
140 MB of memory. The SHA-256 of the uncompressed certificate JSON is
`0ce1bc6d77306573a2eccad01b75168a5faa6d6a058b91c8dcbbd0d1ce6e979a`.

A floating-point linear-programming solver found the multipliers, which were
then rounded to rationals. The solver program is not part of this repository
and has no part in the proof.

## Check the Lean proof

Toolchain: Lean `v4.33.1`, Mathlib `v4.33.1`, and
[lean-formalizations](https://github.com/mysticflounder/lean-formalizations)
at commit `dd46c17`, as pinned in `lean/lake-manifest.json`.

```sh
cd lean
lake exe cache get
lake build
lake env lean CheckAxioms.lean
```

The final statement is

```lean
theorem Esgk.numDistances_gt_one_third_add_kappa_image_gap {n : ℕ} {p : Config n}
    (hp : Function.Injective p) (hgp : InGeneralPosition (Config.toFinset p))
    (hn : 10000 ≤ n) :
    (1 / 3 + 1 / 40) * (n : ℝ) < (NumDistances p : ℝ)
```

With the Mathlib cache, `lake build` uses about 2000 s of CPU time and up to
16 GB of memory. `CheckAxioms.lean` prints its axioms: `propext`, `Classical.choice` and
`Quot.sound`. The proof uses no `sorry` and no `native_decide`. The 60
certificate cases are checked by `decide +kernel`, with the multipliers
scaled to integers at 2²⁰⁰.

The Lean files are copied without change from the research repository in
which the proof was developed. Their comments cite research notes by file
name (`docs/...`). Those notes are not part of this repository. The paper is
the reference for the mathematics.

## License

- The paper text (`paper/`) is licensed under
  [CC BY 4.0](paper/LICENSE).
- The verifier, the certificate and the Lean sources are licensed under the
  [Apache License 2.0](LICENSE).

## Citation

Adam McKenna, *Distinct distances in general position: an explicit
improvement of the bound n/3*, draft, October 2026,
https://github.com/mysticflounder/distinct-distances-general-position.
