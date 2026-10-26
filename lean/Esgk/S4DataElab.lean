/-
Copyright (c) 2026 Adam McKenna. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Adam McKenna
-/

import Esgk.S4RowBuilder

/-!
# Data terms for the `1/40` certificate

The certificate data of Theorem S4 (grid points, bands, chord conditions and the weighted row
tags of each case) is written in the generated Lean files as string literals. The term
elaborators of this file turn such a string into a term:

* `s4_qtree% "n d …"`: an `LTree ℚ` of the rationals `n/d` (one pair per line);
* `s4_ntree3% "a b c …"`: an `LTree (Nat × Nat × Nat)` (one triple per line);
* `s4_ctree% "n₁ d₁ n₂ d₂ n₃ d₃ …"`: an `LTree (ℚ × ℚ × ℚ)` (one triple per line);
* `s4_rows% "z tag args …"`: a balanced `WTree Tag` of weighted row tags (one per line).

The elaborators only build terms; the kernel checks every use of the data. `LTree` stores
the entry of index `i` at the path of the binary digits of `i`, lowest digit first, so that
`LTree.get` is a short kernel computation.
-/

namespace Esgk.S4

open Lean Elab Term Meta

/-- A binary tree with values at the leaves. -/
inductive LTree (α : Type) where
  | leaf (a : α)
  | node (l r : LTree α)

/-- The entry of index `k`: the lowest binary digit of `k` picks the subtree. -/
noncomputable def LTree.get {α : Type} (t : LTree α) : Nat → α :=
  LTree.rec (motive := fun _ ↦ Nat → α) (fun a _ ↦ a)
    (fun _ _ fl fr k ↦ bif Nat.beq (k % 2) 0 then fl (k / 2) else fr (k / 2)) t

private def natE (n : Nat) : Expr := mkRawNatLit n

private def ratE (n d : Nat) : Expr :=
  mkApp2 (Lean.mkConst ``mkRat) (mkApp (Lean.mkConst ``Int.ofNat) (natE n)) (natE d)

private def wordsOf (line : String) : List String :=
  (line.splitOn " ").filter (· ≠ "")

private def linesOf (s : String) : List (List String) :=
  ((s.splitOn "\n").map wordsOf).filter (· ≠ [])

private def nat! (w : String) : Nat := w.toNat!

/-- The tree with `xs[i]` at the path of the binary digits of `i` (lowest first); `xs` is
padded to a power of two with its last entry. -/
private partial def lowBitTree (ty : Expr) (xs : Array Expr) : Expr :=
  let leafE (a : Expr) : Expr := mkApp2 (Lean.mkConst ``LTree.leaf) ty a
  let rec go (ys : Array Expr) : Expr :=
    if ys.size ≤ 1 then leafE ys[0]!
    else
      let ev := (Array.range (ys.size / 2)).map (fun i ↦ ys[2 * i]!)
      let od := (Array.range (ys.size / 2)).map (fun i ↦ ys[2 * i + 1]!)
      mkApp3 (Lean.mkConst ``LTree.node) ty (go ev) (go od)
  Id.run do
    let mut ys := xs
    let mut sz := 1
    while sz < ys.size do sz := 2 * sz
    while ys.size < sz do ys := ys.push ys[ys.size - 1]!
    return go ys

private def prodE (α β : Expr) (a b : Expr) : Expr :=
  mkApp4 (Lean.mkConst ``Prod.mk [Level.zero, Level.zero]) α β a b

private def prodT (α β : Expr) : Expr := mkApp2 (Lean.mkConst ``Prod [Level.zero, Level.zero]) α β

elab "s4_qtree% " s:str : term => do
  let xs := (linesOf s.getString).toArray.map fun ws ↦ ratE (nat! ws[0]!) (nat! ws[1]!)
  return lowBitTree (Lean.mkConst ``Rat) xs

elab "s4_ntree3% " s:str : term => do
  let nat := Lean.mkConst ``Nat
  let xs := (linesOf s.getString).toArray.map fun ws ↦
    prodE nat (prodT nat nat) (natE (nat! ws[0]!))
      (prodE nat nat (natE (nat! ws[1]!)) (natE (nat! ws[2]!)))
  return lowBitTree (prodT nat (prodT nat nat)) xs

elab "s4_ctree% " s:str : term => do
  let q := Lean.mkConst ``Rat
  let xs := (linesOf s.getString).toArray.map fun ws ↦
    prodE q (prodT q q) (ratE (nat! ws[0]!) (nat! ws[1]!))
      (prodE q q (ratE (nat! ws[2]!) (nat! ws[3]!)) (ratE (nat! ws[4]!) (nat! ws[5]!)))
  return lowBitTree (prodT q (prodT q q)) xs

/-- The constructor and the number of arguments of a tag name. -/
private def tagCtor : String → Option (Name × Nat)
  | "amono" => some (``Tag.amono, 2)
  | "xmono" => some (``Tag.xmono, 2)
  | "xpos" => some (``Tag.xpos, 1)
  | "xleq" => some (``Tag.xleq, 1)
  | "a2lo" => some (``Tag.a2lo, 1)
  | "mlo" => some (``Tag.mlo, 0)
  | "layer" => some (``Tag.layer, 2)
  | "pa" => some (``Tag.pa, 4)
  | "pb" => some (``Tag.pb, 4)
  | "pc" => some (``Tag.pc, 4)
  | "aa" => some (``Tag.aa, 5)
  | "jj" => some (``Tag.jj, 5)
  | "ab1" => some (``Tag.ab1, 2)
  | "ab2" => some (``Tag.ab2, 1)
  | "ab3" => some (``Tag.ab3, 3)
  | "ab4" => some (``Tag.ab4, 3)
  | "ab6u" => some (``Tag.ab6u, 2)
  | "ab6v" => some (``Tag.ab6v, 3)
  | "caselo" => some (``Tag.caselo, 0)
  | "casehi" => some (``Tag.casehi, 0)
  | _ => none

private partial def balanced (node : Expr → Expr → Expr) (xs : Array Expr) (lo hi : Nat) :
    Expr :=
  if hi - lo ≤ 1 then xs[lo]!
  else
    let m := (lo + hi) / 2
    node (balanced node xs lo m) (balanced node xs m hi)

elab "s4_rows% " s:str : term => do
  let tagT := Lean.mkConst ``Tag
  let mut leaves : Array Expr := #[]
  for ws in linesOf s.getString do
    let some (c, k) := tagCtor ws[1]! | throwError "unknown tag {ws[1]!}"
    unless ws.length == k + 2 do throwError "tag {ws[1]!}: expected {k} arguments"
    let tag := mkAppN (Lean.mkConst c) ((ws.drop 2).toArray.map fun w ↦ natE (nat! w))
    leaves := leaves.push
      (mkApp3 (Lean.mkConst ``FarkasPackedCheck.WTree.leaf) tagT (natE (nat! ws[0]!)) tag)
  if leaves.isEmpty then throwError "no rows"
  return balanced (mkApp3 (Lean.mkConst ``FarkasPackedCheck.WTree.node) tagT) leaves 0 leaves.size

end Esgk.S4
