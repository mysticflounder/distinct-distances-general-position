/-
Copyright (c) 2026 Adam McKenna. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Adam McKenna
-/

import Esgk.BandTriples
import Esgk.FarkasPackedCheck
import Esgk.PalettePairCountCut
import Esgk.S4RowFacts
import Esgk.ThinAnnulus

/-!
# The rows of the `1/40` certificate

Theorem S4 of `docs/results/explicit-constant-gain-1-40-2026-10-21.md` is proved by an
exact-rational certificate: in each of `60` cases for the mean radial excess `ȳ`, a
nonnegative combination of valid linear rows contradicts the box of the normalized
variables. This file fixes the interface between the certificate data and the geometry:

* `Params`: the certificate parameters (the grid, the bands of Lemma AB and the proved chord
  conditions), as computable functions of an index;
* `Setting`: a configuration with the hypotheses of the contradiction argument
  (`n ≥ 10000` and `D ≤ (1/3 + 1/40) n`), with the derived annulus facts;
* `val`: the normalized variables `𝒜(g) = P(g²ρ²)/n²`, `ℛ(g) = R(g²ρ²)/n`, `D/n`, `1/n` and
  the band variables `W, U, V` of Lemma AB, as a function of a variable index;
* `Tag` and `build`: a row tag and the computable builder that checks the side conditions
  of its row family in `ℚ` and returns the rational row (`QRow`), or `none`;
* `QRow.toRow`: the scaling of a rational row by the common denominator of its entries to
  the integer row of `FarkasPackedCheck`, with its soundness `QRow.toRow_sound`;
* `box`: the box of the normalized variables.

The soundness of each row family (`build … = some r → r.Holds (val P S)`) is proved in
separate files.

Variables: `𝒜(gᵢ) ↦ 2i`, `ℛ(gᵢ) ↦ 2i + 1` (`i < ng`), `D/n ↦ 2 ng`, `1/n ↦ 2 ng + 1`, and
for band `b < nb`: `W_b ↦ 2 ng + 2 + 3b`, `U_b ↦ 2 ng + 3 + 3b`, `V_b ↦ 2 ng + 4 + 3b`.
-/

namespace Esgk.S4

open Finset EuclideanGeometry Esgk.FarkasPackedCheck

/-! ### Constants of the certificate -/

/-- The annulus parameter `s̄ = 41821/50000` (`s̄² ≤ 1 - 12/40 - 4/10⁴`). -/
def sig : ℚ := 41821 / 50000

/-- The annulus width `W = 1 - s̄²`. -/
def wid : ℚ := 1 - sig * sig

/-- The upper bound `1/3 + 1/40` of `D/n` under the contradiction hypothesis. -/
def dbar : ℚ := 1 / 3 + 1 / 40

/-- The upper bound `1/10000` of `1/n`. -/
def eps : ℚ := 1 / 10000

/-! ### Certificate parameters and variables -/

/-- The certificate parameters. -/
structure Params where
  /-- the number of grid points -/
  ng : Nat
  /-- the grid point of an index -/
  g : Nat → ℚ
  /-- the number of bands of Lemma AB -/
  nb : Nat
  /-- the grid index of the threshold `λ` of a band -/
  bLam : Nat → Nat
  /-- the grid index of the lower edge of a band -/
  bLo : Nat → Nat
  /-- the grid index of the upper edge of a band -/
  bHi : Nat → Nat
  /-- the number of proved chord conditions -/
  nch : Nat
  /-- the proved chord conditions `(lo, hi, T)` -/
  ch : Nat → ℚ × ℚ × ℚ

/-- The variable `𝒜(gᵢ)`. -/
def aVar (i : Nat) : Nat := 2 * i

/-- The variable `ℛ(gᵢ)`. -/
def rVar (i : Nat) : Nat := 2 * i + 1

namespace Params

variable (P : Params)

/-- The variable `D/n`. -/
def dV : Nat := 2 * P.ng

/-- The variable `1/n`. -/
def eV : Nat := 2 * P.ng + 1

/-- The variable `W_b` (the number of triples of band `b`, over `n²`). -/
def wV (b : Nat) : Nat := 2 * P.ng + 2 + 3 * b

/-- The variable `U_b` (the base excess of band `b`, over `n²`). -/
def uV (b : Nat) : Nat := 2 * P.ng + 3 + 3 * b

/-- The variable `V_b` (the apex excess of band `b`, over `n²`). -/
def vV (b : Nat) : Nat := 2 * P.ng + 4 + 3 * b

/-- The number of variables. -/
def nv : Nat := 2 * P.ng + 2 + 3 * P.nb

/-- The squared grid point `gᵢ²`. -/
def sq (i : Nat) : ℚ := P.g i * P.g i

/-- The index is a grid index. -/
def ok (i : Nat) : Bool := Nat.blt i P.ng

/-- The index is a band index. -/
def okB (b : Nat) : Bool := Nat.blt b P.nb

end Params

/-! ### The configuration -/

/-- A configuration with the hypotheses of the contradiction argument for Theorem S4. -/
structure Setting (n : ℕ) where
  /-- the points -/
  p : Config n
  hp : Function.Injective p
  hgp : InGeneralPosition (Config.toFinset p)
  hn : 10000 ≤ n
  hD : (NumDistances p : ℝ) ≤ (1 / 3 + 1 / 40) * (n : ℝ)

namespace Setting

variable {n : ℕ} (S : Setting n)

/-- The annulus radius `ρ`. -/
noncomputable def ρ : ℝ := distanceAnnulusRadius S.p

/-- The centroid. -/
noncomputable def c : Point := configCentroid S.p

/-- `P(t)`: the number of ordered pairs of distinct points at squared distance at most `t`. -/
noncomputable def pairs (t : ℝ) : ℝ :=
  ∑ j, ((((univ : Finset (Fin n)).filter (fun q ↦ dist (S.p j) (S.p q) ^ 2 ≤ t)).card : ℝ) - 1)

/-- `R(t)`: the number of palette values `r` with `r² ≤ t`. -/
noncomputable def pal (t : ℝ) : ℝ :=
  (((OrderedDistanceValues S.p).filter (fun r ↦ r ^ 2 ≤ t)).card : ℝ)

/-- The number `D` of distances. -/
noncomputable def D : ℝ := (NumDistancesOrdered S.p : ℝ)

/-- The radial excess `yᵢ = |pᵢ - c|²/ρ² - s̄²`. -/
noncomputable def y (i : Fin n) : ℝ := radialExcess S.p S.c S.ρ (sig : ℝ) i

/-- The mean radial excess `ȳ`. -/
noncomputable def ybar : ℝ := (∑ i, S.y i) / n

/-- The triples of band `b` of Lemma AB. -/
noncomputable def bandTri (P : Params) (b : Nat) : Finset (Fin n × Fin n × Fin n) :=
  apexBandTriples S.p S.ρ ((P.g (P.bLam b) : ℝ) ^ 2) ((P.g (P.bLo b) : ℝ) ^ 2)
    ((P.g (P.bHi b) : ℝ) ^ 2)

/-- The band variables `W_b` (`r = 0`), `U_b` (`r = 1`) and `V_b` (`r = 2`), times `n²`. -/
noncomputable def bandCount (P : Params) (b r : Nat) : ℝ :=
  if r = 0 then ((S.bandTri P b).card : ℝ)
  else if r = 1 then ∑ x ∈ S.bandTri P b, (S.y x.2.1 + S.y x.2.2)
  else ∑ x ∈ S.bandTri P b, S.y x.1

/-- The chord condition `(lo, hi, T)`: every saturated incidence of radius `r` with
`lo²ρ² < r² ≤ hi²ρ²` has two of its points at squared distance at most `Tρ²`. -/
def ChordHolds (lo hi T : ℚ) : Prop :=
  ∀ ir ∈ SaturatedDistanceIncidences S.p, (lo : ℝ) ^ 2 * S.ρ ^ 2 < ir.2 ^ 2 →
    ir.2 ^ 2 ≤ (hi : ℝ) ^ 2 * S.ρ ^ 2 → ∃ a b,
      a ∈ DistanceOutNeighborsAtValue S.p ir.1 ir.2 ∧
      b ∈ DistanceOutNeighborsAtValue S.p ir.1 ir.2 ∧ a < b ∧
      dist (S.p a) (S.p b) ^ 2 ≤ (T : ℝ) * S.ρ ^ 2

theorem n_pos (S : Setting n) : 0 < n := by have := S.hn; omega

theorem hnR (S : Setting n) : (10000 : ℝ) ≤ n := by exact_mod_cast S.hn

theorem hcap : ∀ i r, (DistanceOutNeighborsAtValue S.p i r).card ≤ 3 :=
  distanceOutNeighborsAtValue_card_le_three_of_noFourCocircular S.hp S.hgp.2

theorem hρ : 0 < S.ρ := distanceAnnulusRadius_pos S.hp (by have := S.hn; omega) S.hcap

theorem hρsq : S.ρ ^ 2 = distanceAnnulusRadiusSq S.p :=
  (distanceAnnulusRadius_bounds S.n_pos S.hcap).2.1

theorem hD_ordered : S.D ≤ (1 / 3 + 1 / 40) * (n : ℝ) := by
  have h := S.hD
  rw [numDistances_eq_numDistancesOrdered_of_injective S.hp] at h
  exact h

theorem n_le_three_D : (n : ℝ) ≤ 3 * S.D :=
  numPoints_le_three_mul_numDistancesOrdered_real S.hp S.hgp (by have := S.hn; omega)

theorem hhi (i : Fin n) : dist (S.p i) S.c ^ 2 ≤ S.ρ ^ 2 := by
  rw [S.hρsq]
  exact (thin_annulus_sq_of_noFourCocircular S.hp S.hgp.2 S.n_pos S.hD i).2

theorem sig_cast : ((sig : ℚ) : ℝ) = 41821 / 50000 := by
  norm_num [sig]

theorem hlo (i : Fin n) : (sig : ℝ) ^ 2 * S.ρ ^ 2 ≤ dist (S.p i) S.c ^ 2 := by
  have h := (thin_annulus_sq_of_noFourCocircular S.hp S.hgp.2 S.n_pos S.hD i).1
  rw [← S.hρsq] at h
  have h4n : 4 / (n : ℝ) ≤ 4 / 10000 :=
    div_le_div_of_nonneg_left (by norm_num) (by norm_num) S.hnR
  have hσ : (sig : ℝ) ^ 2 ≤ 1 - 12 * (1 / 40 : ℝ) - 4 / (n : ℝ) := by
    rw [sig_cast]
    linarith
  have hρ2 : 0 ≤ S.ρ ^ 2 := sq_nonneg _
  have := mul_le_mul_of_nonneg_right hσ hρ2
  show (sig : ℝ) ^ 2 * S.ρ ^ 2 ≤ dist (S.p i) (configCentroid S.p) ^ 2
  linarith

theorem hsum : ∑ i, (S.p i - S.c) = 0 := sum_sub_configCentroid S.p S.n_pos

end Setting

/-! ### The valuation -/

/-- The normalized variables of the certificate, by variable index. -/
noncomputable def val {n : ℕ} (P : Params) (S : Setting n) (k : Nat) : ℝ :=
  if k < 2 * P.ng then
    if k % 2 = 0 then S.pairs ((P.g (k / 2) : ℝ) ^ 2 * S.ρ ^ 2) / ((n : ℝ) * n)
    else S.pal ((P.g (k / 2) : ℝ) ^ 2 * S.ρ ^ 2) / n
  else if k = 2 * P.ng then S.D / n
  else if k = 2 * P.ng + 1 then 1 / n
  else S.bandCount P ((k - (2 * P.ng + 2)) / 3) ((k - (2 * P.ng + 2)) % 3) / ((n : ℝ) * n)

section Valuation

variable {n : ℕ} (P : Params) (S : Setting n)

theorem val_aVar {i : Nat} (hi : i < P.ng) :
    val P S (aVar i) = S.pairs ((P.g i : ℝ) ^ 2 * S.ρ ^ 2) / ((n : ℝ) * n) := by
  unfold val aVar
  rw [if_pos (by omega), if_pos (by omega), show 2 * i / 2 = i by omega]

theorem val_rVar {i : Nat} (hi : i < P.ng) :
    val P S (rVar i) = S.pal ((P.g i : ℝ) ^ 2 * S.ρ ^ 2) / n := by
  unfold val rVar
  rw [if_pos (by omega), if_neg (by omega), show (2 * i + 1) / 2 = i by omega]

theorem val_dV : val P S P.dV = S.D / n := by
  unfold val Params.dV
  rw [if_neg (by omega), if_pos rfl]

theorem val_eV : val P S P.eV = 1 / n := by
  unfold val Params.eV
  rw [if_neg (by omega), if_neg (by omega), if_pos rfl]

theorem val_wV (b : Nat) :
    val P S (P.wV b) = ((S.bandTri P b).card : ℝ) / ((n : ℝ) * n) := by
  unfold val Params.wV
  rw [if_neg (by omega), if_neg (by omega), if_neg (by omega),
    show (2 * P.ng + 2 + 3 * b - (2 * P.ng + 2)) / 3 = b by omega,
    show (2 * P.ng + 2 + 3 * b - (2 * P.ng + 2)) % 3 = 0 by omega]
  simp [Setting.bandCount]

theorem val_uV (b : Nat) :
    val P S (P.uV b) = (∑ x ∈ S.bandTri P b, (S.y x.2.1 + S.y x.2.2)) / ((n : ℝ) * n) := by
  unfold val Params.uV
  rw [if_neg (by omega), if_neg (by omega), if_neg (by omega),
    show (2 * P.ng + 3 + 3 * b - (2 * P.ng + 2)) / 3 = b by omega,
    show (2 * P.ng + 3 + 3 * b - (2 * P.ng + 2)) % 3 = 1 by omega]
  simp [Setting.bandCount]

theorem val_vV (b : Nat) :
    val P S (P.vV b) = (∑ x ∈ S.bandTri P b, S.y x.1) / ((n : ℝ) * n) := by
  unfold val Params.vV
  rw [if_neg (by omega), if_neg (by omega), if_neg (by omega),
    show (2 * P.ng + 4 + 3 * b - (2 * P.ng + 2)) / 3 = b by omega,
    show (2 * P.ng + 4 + 3 * b - (2 * P.ng + 2)) % 3 = 2 by omega]
  simp [Setting.bandCount]

end Valuation

/-! ### Rational rows -/

/-- A rational row `∑_{(k, q) ∈ ent} q * v k ≤ rhs`. -/
structure QRow where
  /-- the entries `(variable, coefficient)`; an index may repeat -/
  ent : List (Nat × ℚ)
  /-- the right-hand side -/
  rhs : ℚ

/-- The left-hand side of a rational row at `v`. -/
noncomputable def QRow.lhs (r : QRow) (v : Nat → ℝ) : ℝ :=
  (r.ent.map fun e ↦ (e.2 : ℝ) * v e.1).sum

/-- The point `v` satisfies the rational row `r`. -/
def QRow.Holds (r : QRow) (v : Nat → ℝ) : Prop := r.lhs v ≤ r.rhs

/-! ### Scaling a rational row to an integer row -/

/-- The lcm of the denominators of the entries and of the right-hand side. -/
noncomputable def QRow.den (r : QRow) : Nat :=
  List.rec (motive := fun _ ↦ Nat) r.rhs.den (fun (e : Nat × ℚ) _ ih ↦ Nat.lcm e.2.den ih) r.ent

/-- Every entry times `L` is an integer. -/
noncomputable def allInt (L : Nat) (E : List (Nat × ℚ)) : Bool :=
  List.rec (motive := fun _ ↦ Bool) true (fun (e : Nat × ℚ) _ ih ↦ Nat.beq ((e.2 : ℚ) * (L : ℚ)).den 1 && ih) E

/-- The positive entries times `L`. -/
noncomputable def posPart (L : Nat) (E : List (Nat × ℚ)) : List (Nat × Nat) :=
  List.rec (motive := fun _ ↦ List (Nat × Nat)) []
    (fun (e : Nat × ℚ) _ ih ↦
      bif decide (0 < e.2) then (e.1, ((e.2 : ℚ) * (L : ℚ)).num.toNat) :: ih else ih) E

/-- The negative entries times `-L`. -/
noncomputable def negPart (L : Nat) (E : List (Nat × ℚ)) : List (Nat × Nat) :=
  List.rec (motive := fun _ ↦ List (Nat × Nat)) []
    (fun (e : Nat × ℚ) _ ih ↦
      bif decide (e.2 < 0) then (e.1, (-((e.2 : ℚ) * (L : ℚ))).num.toNat) :: ih else ih) E

/-- The integer row `L * r` for the common denominator `L` of `r`, or `none` if the scaling
fails (it does not fail, but the checker does not need that). -/
noncomputable def QRow.toRow (r : QRow) : Option Row :=
  bif Nat.blt 0 r.den && allInt r.den r.ent && Nat.beq (r.rhs * r.den).den 1 then
    some ⟨posPart r.den r.ent, negPart r.den r.ent, (r.rhs * r.den).num⟩
  else none


private theorem allInt_cons (L : Nat) (e : Nat × ℚ) (E : List (Nat × ℚ)) :
    allInt L (e :: E) = (Nat.beq (e.2 * L).den 1 && allInt L E) := rfl

private theorem posPart_cons (L : Nat) (e : Nat × ℚ) (E : List (Nat × ℚ)) :
    posPart L (e :: E) =
      bif decide (0 < e.2) then (e.1, (e.2 * L).num.toNat) :: posPart L E else posPart L E :=
  rfl

private theorem negPart_cons (L : Nat) (e : Nat × ℚ) (E : List (Nat × ℚ)) :
    negPart L (e :: E) =
      bif decide (e.2 < 0) then (e.1, (-(e.2 * L)).num.toNat) :: negPart L E
      else negPart L E :=
  rfl

private theorem entrySum_cons (x : Nat × Nat) (E : List (Nat × Nat)) (v : Nat → ℝ) :
    entrySum (x :: E) v = (x.2 : ℝ) * v x.1 + entrySum E v := by
  simp [entrySum]

private theorem intCast_num_of_den_eq_one {q : ℚ} (hd : q.den = 1) :
    ((q.num : ℤ) : ℝ) = (q : ℝ) := by
  rw [← Rat.cast_intCast, (Rat.den_eq_one_iff q).1 hd]

private theorem natCast_toNat_num {q : ℚ} (hd : q.den = 1) (h0 : 0 ≤ q) :
    ((q.num.toNat : ℕ) : ℝ) = (q : ℝ) := by
  rw [← Int.cast_natCast, Int.toNat_of_nonneg (Rat.num_nonneg.2 h0),
    intCast_num_of_den_eq_one hd]

private theorem entrySum_pos_sub_neg (L : Nat) (v : Nat → ℝ) :
    ∀ E : List (Nat × ℚ), allInt L E = true →
      entrySum (posPart L E) v - entrySum (negPart L E) v =
        (L : ℝ) * (E.map fun e ↦ (e.2 : ℝ) * v e.1).sum
  | [], _ => by simp [entrySum, posPart, negPart]
  | e :: E, h => by
    rw [allInt_cons, Bool.and_eq_true] at h
    have ih := entrySum_pos_sub_neg L v E h.2
    have hden : (e.2 * L).den = 1 := Nat.eq_of_beq_eq_true h.1
    rw [posPart_cons, negPart_cons, List.map_cons, List.sum_cons]
    have hL : (0 : ℚ) ≤ L := Nat.cast_nonneg L
    rcases lt_trichotomy e.2 0 with hneg | hzero | hpos
    · rw [decide_eq_false (show ¬ (0 < e.2) by linarith), decide_eq_true hneg, cond_false,
        cond_true, entrySum_cons]
      have hc := natCast_toNat_num (q := -(e.2 * L)) (by rw [Rat.neg_den]; exact hden)
        (neg_nonneg.2 (mul_nonpos_of_nonpos_of_nonneg hneg.le hL))
      push_cast at hc
      rw [hc]
      linear_combination ih
    · simp only [hzero, lt_irrefl, decide_false, cond_false, Rat.cast_zero, zero_mul, zero_add]
      exact ih
    · rw [decide_eq_true hpos, decide_eq_false (show ¬ (e.2 < 0) by linarith), cond_true,
        cond_false, entrySum_cons]
      have hc := natCast_toNat_num hden (mul_nonneg hpos.le hL)
      push_cast at hc
      rw [hc]
      linear_combination ih

/-- Soundness of the scaling: the integer row holds wherever the rational row holds. -/
theorem QRow.toRow_sound {r : QRow} {R : Row} {v : Nat → ℝ} (h : r.toRow = some R)
    (hr : r.Holds v) : R.Holds v := by
  unfold QRow.toRow at h
  cases hc : (Nat.blt 0 r.den && allInt r.den r.ent && Nat.beq (r.rhs * r.den).den 1) with
  | false => rw [hc, cond_false] at h; exact absurd h (by simp)
  | true =>
    rw [hc, cond_true] at h
    simp only [Bool.and_eq_true] at hc
    obtain ⟨⟨hL, hint⟩, hrhs⟩ := hc
    injection h with h
    subst h
    have hL' : (0 : ℝ) < r.den := by exact_mod_cast (by simpa [Nat.blt_eq] using hL : 0 < r.den)
    unfold Row.Holds Row.lhs
    rw [entrySum_pos_sub_neg _ _ _ hint, intCast_num_of_den_eq_one (Nat.eq_of_beq_eq_true hrhs)]
    push_cast
    unfold QRow.Holds QRow.lhs at hr
    have := mul_le_mul_of_nonneg_left hr hL'.le
    linarith

/-! ### Row tags -/

/-- A row tag: the row family and its grid indices, band indices and rational witnesses
(`(num, den)` pairs). -/
inductive Tag where
  /-- `𝒜(gᵢ) ≤ 𝒜(gⱼ)` -/
  | amono (i j : Nat)
  /-- `3ℛ - 𝒜` is nondecreasing from `gᵢ` to `gⱼ` -/
  | xmono (i j : Nat)
  /-- `𝒜(gᵢ) ≤ 3ℛ(gᵢ)` -/
  | xpos (i : Nat)
  /-- `3ℛ(gᵢ) - 𝔮 ≤ 𝒜(gᵢ)` -/
  | xleq (i : Nat)
  /-- `1 - 1/n ≤ 𝒜(gᵢ)` for `gᵢ ≥ 2` -/
  | a2lo (i : Nat)
  /-- Lemma M (lower Riemann sum over the whole grid) -/
  | mlo
  /-- Lemma K over the grid points `gᵢ, …, gⱼ` -/
  | layer (i j : Nat)
  /-- Lemma P″(a) at `ρ = g_iρ`, `r′ = g_ir`, with `g_ih = ρ/2` and chord evidence `e` -/
  | pa (iρ ir ih e : Nat)
  /-- Lemma P″(b) -/
  | pb (iρ ir ih e : Nat)
  /-- Lemma P″(c) -/
  | pc (iρ ir ih e : Nat)
  /-- Lemma A″ at `λ = g_il`, `w = g_iw`, `h = g_ih`, with the witness `r₄ ≥ √(4 - λ²)` -/
  | aa (il iw ih r4n r4d : Nat)
  /-- Lemma J at `λ = g_il`, `g₂ = g_ig2`, `g_σ = g_igs`, with the witness
  `u ≤ √(s̄² - h₀²)` -/
  | jj (il ig2 igs un ud : Nat)
  /-- Lemma AB, row AB1, over the bands `b₀, …, b₀ + m - 1` of one threshold -/
  | ab1 (b0 m : Nat)
  /-- Lemma AB, row AB2, on band `b` -/
  | ab2 (b : Nat)
  /-- Lemma AB, row AB3, on band `b` at `t = tn/td` -/
  | ab3 (b tn td : Nat)
  /-- Lemma AB, row AB4, on band `b` at `t = tn/td` -/
  | ab4 (b tn td : Nat)
  /-- Lemma AB, row AB6u, over the bands `b₀, …, b₀ + m - 1` -/
  | ab6u (b0 m : Nat)
  /-- Lemma AB, row AB6v, over the bands `b₀, …, b₀ + m - 1`, with `h = g_ih` -/
  | ab6v (b0 m ih : Nat)
  /-- the upper case row (`Y_lo ≤ ȳ` with the upper Riemann sum) -/
  | caselo
  /-- the lower case row (`ȳ ≤ Y_hi` with the lower Riemann sum) -/
  | casehi

/-! ### The row builder -/

namespace Params

variable (P : Params)

/-- `gᵢ² ≤ gᵢ₊₁²` for the `c` cells starting at `k`. -/
def gridMono : Nat → Nat → Bool
  | 0, _ => true
  | c + 1, k => decide (P.sq k ≤ P.sq (k + 1)) && gridMono c (k + 1)

/-- The entries `f k, f (k + 1), …, f (k + c - 1)`, concatenated. -/
def cells (f : Nat → List (Nat × ℚ)) : Nat → Nat → List (Nat × ℚ)
  | 0, _ => []
  | c + 1, k => f k ++ cells f c (k + 1)

/-- The bands `b₀, …, b₀ + m - 1` share the threshold of `b₀`, have nondecreasing edges, and
the upper edge of each band is the lower edge of the next. -/
def chain (b0 : Nat) : Nat → Bool
  | 0 => true
  | m + 1 => chain b0 m && Nat.beq (P.bLam (b0 + m)) (P.bLam b0) &&
      decide (P.sq (P.bLo (b0 + m)) ≤ P.sq (P.bHi (b0 + m))) &&
      (Nat.beq m 0 || Nat.beq (P.bHi (b0 + m - 1)) (P.bLo (b0 + m)))

/-- The chord evidence `e` for the condition `(lo, hi, T)`: `0` is Step L0, `1` is the radius
bound, `k + 2` is the proved chord condition `k`. -/
def chordOK (lo hi T : ℚ) : Nat → Bool
  | 0 => decide ((1 - sig * sig) * (1 + sig) ≤ sig * T)
  | 1 => decide ((1 - sig * sig) * (1 - sig * sig) * (hi * hi) ≤ sig * sig * (T * T)) &&
      decide (0 ≤ T)
  | k + 2 => Nat.blt k P.nch && decide (0 ≤ (P.ch k).1) && decide ((P.ch k).1 ≤ lo) &&
      decide (0 ≤ hi) && decide (hi ≤ (P.ch k).2.1) && decide ((P.ch k).2.2 ≤ T)

/-- Step L0 at the squared threshold `l`. -/
def stepL0 (l : ℚ) : Bool := decide ((1 - sig * sig) * (sig + 1) ≤ sig * l)

/-- The rows `𝒜(gᵢ) ≤ 𝒜(gⱼ)`. -/
def rowAmono (i j : Nat) : Option QRow :=
  if P.ok i && P.ok j && decide (P.sq i ≤ P.sq j) then
    some ⟨[(aVar i, 1), (aVar j, -1)], 0⟩
  else none

/-- The rows `3ℛ(gᵢ) - 𝒜(gᵢ) ≤ 3ℛ(gⱼ) - 𝒜(gⱼ)`. -/
def rowXmono (i j : Nat) : Option QRow :=
  if P.ok i && P.ok j && decide (P.sq i ≤ P.sq j) then
    some ⟨[(rVar i, 3), (aVar i, -1), (rVar j, -3), (aVar j, 1)], 0⟩
  else none

/-- The rows `𝒜(gᵢ) ≤ 3ℛ(gᵢ)`. -/
def rowXpos (i : Nat) : Option QRow :=
  if P.ok i then some ⟨[(aVar i, 1), (rVar i, -3)], 0⟩ else none

/-- The rows `3ℛ(gᵢ) - 𝒜(gᵢ) - 3 D/n - 1/n ≤ -1`. -/
def rowXleq (i : Nat) : Option QRow :=
  if P.ok i then some ⟨[(rVar i, 3), (aVar i, -1), (P.dV, -3), (P.eV, -1)], -1⟩ else none

/-- The row `-𝒜(gᵢ) - 1/n ≤ -1` for `gᵢ² ≥ 4`. -/
def rowA2lo (i : Nat) : Option QRow :=
  if P.ok i && decide (4 ≤ P.sq i) then some ⟨[(aVar i, -1), (P.eV, -1)], -1⟩ else none

/-- The entries of one cell of Lemma M. -/
def mloCell (k : Nat) : List (Nat × ℚ) :=
  [(aVar k, -((P.sq (k + 1) - P.sq k) / 2)), (rVar k, 6 * ((P.sq (k + 1) - P.sq k) / 2))]

/-- The row of Lemma M: `∑ w (6ℛ(a) - 𝒜(a)) - 12 D/n - 2/n ≤ -3`. -/
def rowMlo : Option QRow :=
  if Nat.blt 0 P.ng && decide (0 ≤ P.sq 0) && P.gridMono (P.ng - 1) 0 &&
      decide (P.sq (P.ng - 1) ≤ 4) then
    some ⟨cells P.mloCell (P.ng - 1) 0 ++ [(P.dV, -12), (P.eV, -2)], -3⟩
  else none

/-- The entries of one cell of Lemma K. -/
def layerCell (k : Nat) : List (Nat × ℚ) := [(aVar k, P.sq (k + 1) - P.sq k)]

/-- The row of Lemma K over `gᵢ, …, gⱼ`: `∑ (b² - a²) 𝒜(a) ≤ 2s̄`. -/
def rowLayer (i j : Nat) : Option QRow :=
  if P.ok j && Nat.ble i j && P.gridMono (j - i) i &&
      decide (P.sq j ≤ 2 * sig * sig + 2 * sig) then
    some ⟨cells P.layerCell (j - i) i, 2 * sig⟩
  else none

/-- The side conditions of Lemma P″ at `ρ = g_iρ`, `r′ = g_ir`, `ρ/2 = g_ih`. -/
def palOK (iρ ir ih e : Nat) : Bool :=
  P.ok iρ && P.ok ir && P.ok ih && decide (0 ≤ P.g ih) && decide (4 * P.sq ih ≤ P.sq iρ) &&
    decide (P.sq ih ≤ P.sq ir) && P.chordOK (P.g ih) (P.g ir) (P.sq iρ) e

/-- The common part of the rows of Lemma P″: `𝒜(r′) - 2ℛ(r′) - 𝒜(ρ) ≤ 0`, or for `r′ ≥ 2`,
`-1/n - 2 D/n - 𝒜(ρ) ≤ -1`. -/
def palBase (iρ ir : Nat) : List (Nat × ℚ) × ℚ :=
  bif decide (4 ≤ P.sq ir) then ([(P.eV, -1), (P.dV, -2), (aVar iρ, -1)], -1)
  else ([(aVar ir, 1), (rVar ir, -2), (aVar iρ, -1)], 0)

/-- The row of Lemma P″(a). -/
def rowPa (iρ ir ih e : Nat) : Option QRow :=
  if P.palOK iρ ir ih e then some ⟨(P.palBase iρ ir).1, (P.palBase iρ ir).2⟩ else none

/-- The row of Lemma P″(b): Lemma P″(a) plus `𝒜(ρ/2) - ℛ(ρ/2)`. -/
def rowPb (iρ ir ih e : Nat) : Option QRow :=
  if P.palOK iρ ir ih e then
    some ⟨(P.palBase iρ ir).1 ++ [(aVar ih, 1), (rVar ih, -1)], (P.palBase iρ ir).2⟩
  else none

/-- The row of Lemma P″(c): Lemma P″(a) minus `𝒜(ρ/2) - 2ℛ(ρ/2)`. -/
def rowPc (iρ ir ih e : Nat) : Option QRow :=
  if P.palOK iρ ir ih e then
    some ⟨(P.palBase iρ ir).1 ++ [(aVar ih, -1), (rVar ih, 2)], (P.palBase iρ ir).2⟩
  else none

/-- The row of Lemma A″: `-2𝒜(λ) - (4/3)𝒜(h) + (4/3)𝒜(w) - 6 D/n - 4/n ≤ -4`, with
`r₄ = r4n/r4d`. -/
def rowAA (il iw ih r4n r4d : Nat) : Option QRow :=
  let r4 := mkRat r4n r4d
  if P.ok il && P.ok iw && P.ok ih && decide (0 ≤ r4) && decide (4 - P.sq il ≤ r4 * r4) &&
      decide (4 - P.sq il ≤ 4 * (sig * sig)) && stepL0 (P.sq il) &&
      decide (0 ≤ P.sq iw) && decide (P.sq iw ≤ sig * sig + 1 - sig * r4) &&
      decide (2 + r4 ≤ P.sq ih) then
    some ⟨[(aVar il, -2), (aVar ih, -4 / 3), (aVar iw, 4 / 3), (P.dV, -6), (P.eV, -4)], -4⟩
  else none

/-- The row of Lemma J: `𝒜(g₂) - 𝒜(λ) - 𝒜(g_σ) - 6 D/n - 2/n ≤ -2`, with `h₀ = W/(2λ)`
(`H = h₀²`) and `u = un/ud`. -/
def rowJ (il ig2 igs un ud : Nat) : Option QRow :=
  let u := mkRat un ud
  let H := wid * wid / (4 * P.sq il)
  if P.ok il && P.ok ig2 && P.ok igs && decide (0 < P.sq il) &&
      decide (wid * wid ≤ 4 * P.sq il * H) && decide (P.sq ig2 < 4 * (sig * sig - H)) &&
      decide (0 ≤ u) && decide (u * u ≤ sig * sig - H) && decide (0 ≤ P.g igs) &&
      decide (1 - sig * sig ≤ P.g igs * P.g igs + 2 * P.g igs * u) then
    some ⟨[(aVar ig2, 1), (aVar il, -1), (aVar igs, -1), (P.dV, -6), (P.eV, -2)], -2⟩
  else none

/-- The entries `(sel (b₀ + j), q)` for `j < m`. -/
def bandEnt (sel : Nat → Nat) (q : ℚ) (b0 : Nat) : Nat → List (Nat × ℚ)
  | 0 => []
  | m + 1 => bandEnt sel q b0 m ++ [(sel (b0 + m), q)]

/-- The row AB1: `-2𝒜(λ) - 6 D/n - 4/n - ∑ⱼ W_j ≤ -4`. -/
def rowAB1 (b0 m : Nat) : Option QRow :=
  if Nat.blt 0 m && P.okB (b0 + m - 1) && P.chain b0 m && P.ok (P.bLam b0) &&
      decide (4 * P.sq (P.bLo b0) ≤ P.sq (P.bLam b0)) &&
      decide (4 ≤ P.sq (P.bHi (b0 + m - 1))) then
    some ⟨[(aVar (P.bLam b0), -2), (P.dV, -6), (P.eV, -4)] ++ bandEnt P.wV (-1) b0 m, -4⟩
  else none

/-- The row AB2: `W_b ≤ (4/3)(𝒜(e_hi) - 𝒜(e_lo))`. -/
def rowAB2 (b : Nat) : Option QRow :=
  if P.okB b && P.ok (P.bLo b) && P.ok (P.bHi b) && stepL0 (P.sq (P.bLam b)) &&
      decide (P.sq (P.bLo b) ≤ P.sq (P.bHi b)) then
    some ⟨[(P.wV b, 1), (aVar (P.bHi b), -4 / 3), (aVar (P.bLo b), 4 / 3)], 0⟩
  else none

/-- The row AB3 at `t`: `(2e_lo² - ((4 + t)s̄² + (4s̄² - λ²)/t)) W_b - (1 + 2/t) U_b
- (2 + t) V_b ≤ 0`. -/
def rowAB3 (b tn td : Nat) : Option QRow :=
  let t := mkRat tn td
  if P.okB b && decide (0 < t) then
    some ⟨[(P.wV b, 2 * P.sq (P.bLo b) -
        ((4 + t) * (sig * sig) + (4 * (sig * sig) - P.sq (P.bLam b)) / t)),
      (P.uV b, -(1 + 2 / t)), (P.vV b, -(2 + t))], 0⟩
  else none

/-- `max x 0`. -/
def posPartQ (x : ℚ) : ℚ := bif decide (x ≤ 0) then 0 else x

/-- The row AB4 at `t`: `((4 - t)s̄² - (4s̄² - λ²)/t - 2e_hi²) W_b - max(2/t - 1, 0) U_b
- max(t - 2, 0) V_b ≤ 0`. -/
def rowAB4 (b tn td : Nat) : Option QRow :=
  let t := mkRat tn td
  if P.okB b && decide (0 < t) then
    some ⟨[(P.wV b, ((4 - t) * (sig * sig) - (4 * (sig * sig) - P.sq (P.bLam b)) / t) -
        2 * P.sq (P.bHi b)),
      (P.uV b, -posPartQ (2 / t - 1)), (P.vV b, -posPartQ (t - 2))], 0⟩
  else none

/-- The row AB6u for the case `ȳ ≤ Y_hi`: `∑ⱼ U_j - 12 Y_hi D/n + 12 Y_hi ℛ(λ) ≤ 0`. -/
def rowAB6u (Yhi : ℚ) (b0 m : Nat) : Option QRow :=
  if Nat.blt 0 m && P.okB (b0 + m - 1) && P.chain b0 m && P.ok (P.bLam b0) then
    some ⟨bandEnt P.uV 1 b0 m ++ [(P.dV, -12 * Yhi), (rVar (P.bLam b0), 12 * Yhi)], 0⟩
  else none

/-- The row AB6v for the case `ȳ ≤ Y_hi`: `∑ⱼ V_j - 4 Y_hi D/n + 4 Y_hi ℛ(h) ≤ 0`. -/
def rowAB6v (Yhi : ℚ) (b0 m ih : Nat) : Option QRow :=
  if Nat.blt 0 m && P.okB (b0 + m - 1) && P.chain b0 m && P.ok ih &&
      stepL0 (P.sq (P.bLam b0)) && decide (4 * P.sq ih ≤ P.sq (P.bLam b0)) then
    some ⟨bandEnt P.vV 1 b0 m ++ [(P.dV, -4 * Yhi), (rVar ih, 4 * Yhi)], 0⟩
  else none

/-- The entries of one cell of the upper case row. -/
def upCell (k : Nat) : List (Nat × ℚ) :=
  [(rVar (k + 1), -6 * ((P.sq (k + 1) - P.sq k) / 2)),
    (aVar (k + 1), 2 * ((P.sq (k + 1) - P.sq k) / 2))]

/-- The upper case row for `Y_lo ≤ ȳ`: `12 D/n + 4/n - 2 Up ≤ W - Y_lo + 4`. -/
def rowCaseLo (Ylo : ℚ) : Option QRow :=
  if Nat.blt 0 P.ng && decide (0 ≤ P.sq 0) && P.gridMono (P.ng - 1) 0 &&
      decide (4 ≤ P.sq (P.ng - 1)) then
    some ⟨[(P.dV, 12), (P.eV, 4)] ++ cells P.upCell (P.ng - 1) 0 ++
        [(rVar 0, -6 * (P.sq 0 / 2)), (aVar 0, 2 * (P.sq 0 / 2))], wid - Ylo + 4⟩
  else none

/-- The entries of one cell of the lower case row. -/
def loCell (k : Nat) : List (Nat × ℚ) :=
  [(rVar k, 6 * ((P.sq (k + 1) - P.sq k) / 2)), (aVar k, -2 * ((P.sq (k + 1) - P.sq k) / 2))]

/-- The lower case row for `ȳ ≤ Y_hi`: `-12 D/n - 4/n + 2 Lo ≤ Y_hi - W - 4`. -/
def rowCaseHi (Yhi : ℚ) : Option QRow :=
  if Nat.blt 0 P.ng && decide (0 ≤ P.sq 0) && P.gridMono (P.ng - 1) 0 &&
      decide (P.sq (P.ng - 1) ≤ 4) then
    some ⟨[(P.dV, -12), (P.eV, -4)] ++ cells P.loCell (P.ng - 1) 0, Yhi - wid - 4⟩
  else none

/-- The row of a tag in the case `Y_lo ≤ ȳ ≤ Y_hi`, or `none` if a side condition fails. -/
def build (Ylo Yhi : ℚ) : Tag → Option QRow
  | .amono i j => P.rowAmono i j
  | .xmono i j => P.rowXmono i j
  | .xpos i => P.rowXpos i
  | .xleq i => P.rowXleq i
  | .a2lo i => P.rowA2lo i
  | .mlo => P.rowMlo
  | .layer i j => P.rowLayer i j
  | .pa iρ ir ih e => P.rowPa iρ ir ih e
  | .pb iρ ir ih e => P.rowPb iρ ir ih e
  | .pc iρ ir ih e => P.rowPc iρ ir ih e
  | .aa il iw ih r4n r4d => P.rowAA il iw ih r4n r4d
  | .jj il ig2 igs un ud => P.rowJ il ig2 igs un ud
  | .ab1 b0 m => P.rowAB1 b0 m
  | .ab2 b => P.rowAB2 b
  | .ab3 b tn td => P.rowAB3 b tn td
  | .ab4 b tn td => P.rowAB4 b tn td
  | .ab6u b0 m => P.rowAB6u Yhi b0 m
  | .ab6v b0 m ih => P.rowAB6v Yhi b0 m ih
  | .caselo => P.rowCaseLo Ylo
  | .casehi => P.rowCaseHi Yhi

/-- The integer row of a tag. -/
noncomputable def buildInt (Ylo Yhi : ℚ) (t : Tag) : Option Row :=
  (P.build Ylo Yhi t).bind QRow.toRow

/-- The box of the normalized variables: `𝒜 ∈ [0, 1]`, `ℛ ∈ [0, d̄]`,
`D/n ∈ [1/3 - ε/3, d̄]`, `1/n ∈ [0, ε]`, `W ∈ [0, 2]`, `U, V ∈ [0, 4]`. -/
def box : Box where
  nv := P.nv
  lo k := bif Nat.beq k (2 * P.ng) then 1 / 3 - eps / 3 else 0
  hi k := bif Nat.blt k (2 * P.ng) then (bif Nat.beq (k % 2) 0 then 1 else dbar)
    else bif Nat.beq k (2 * P.ng) then dbar
    else bif Nat.beq k (2 * P.ng + 1) then eps
    else bif Nat.beq ((k - (2 * P.ng + 2)) % 3) 0 then 2 else 4

/-- The proved chord conditions hold for `S`. -/
def ChordsHold {n : ℕ} (S : Setting n) : Prop :=
  ∀ k < P.nch, S.ChordHolds (P.ch k).1 (P.ch k).2.1 (P.ch k).2.2

end Params

end Esgk.S4
