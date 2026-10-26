import Esgk.S4ConstantGain

/-! Run `lake env lean CheckAxioms.lean` after `lake build`. Expected output:
`[propext, Classical.choice, Quot.sound]`. -/

#check @Esgk.numDistances_gt_one_third_add_kappa_image_gap
#print axioms Esgk.numDistances_gt_one_third_add_kappa_image_gap
