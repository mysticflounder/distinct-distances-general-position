"""Refresh data/certificate_cases.json from a verifier run.

Usage:
  extract_data.py certificate <verifier output .txt>
      parses the output of verifier/verify_constant_gain_1_40.py into
      data/certificate_cases.json.

data/profile_1_39.json holds a floating-point LP feasible point at
delta = 1/39 (ybar = 0.154).  The solver that produced it is not part of
this repository, so that file has no refresh step here.
"""
import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.join(os.path.dirname(HERE), "data")


def certificate(src):
    sbar = Fraction(41821, 50000)
    omega = 1 - sbar * sbar
    cases = []
    smallest = None
    verified = False
    with open(src) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("case ["):
                # case [0.1502,0.1552]: rows 1936, gap = 3.2604e-06 ok
                head, tail = line.split("]:", 1)
                lo, hi = head[len("case ["):].split(",")
                fields = tail.replace(",", " ").split()
                rows = int(fields[fields.index("rows") + 1])
                gap = float(fields[fields.index("=") + 1])
                status = fields[-1]
                k = len(cases)
                cases.append({
                    "k": k,
                    "ylo": float(omega * k / 60), "yhi": float(omega * (k + 1) / 60),
                    "ylo_printed": float(lo), "yhi_printed": float(hi),
                    "rows": rows, "gap": gap, "status": status,
                })
            elif line.startswith("smallest gap"):
                smallest = float(line.split("=")[1])
            elif line.startswith("VERIFIED"):
                verified = True
    assert len(cases) == 60, len(cases)
    for c in cases:
        assert abs(c["ylo"] - c["ylo_printed"]) < 6e-5, c
        assert abs(c["yhi"] - c["yhi_printed"]) < 6e-5, c
        assert c["status"] == "ok", c
    out = {
        "source": "verifier/verify_constant_gain_1_40.py "
                  "certificate/constant-gain-1-40-K160-c60.json.gz",
        "omega": float(omega), "verified": verified,
        "smallest_gap": smallest, "total_rows": sum(c["rows"] for c in cases),
        "cases": cases,
    }
    path = os.path.join(DATADIR, "certificate_cases.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote", path, "cases", len(cases), "rows", out["total_rows"],
          "smallest gap", smallest, "verified", verified)


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] != "certificate":
        sys.exit(__doc__)
    certificate(sys.argv[2])
