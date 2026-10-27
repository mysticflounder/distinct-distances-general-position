#!/bin/bash
# Build every figure of the paper as a vector PDF in this directory.
#
# Usage:
#   ./build.sh                  build all figures from src/ and data/
#   ./build.sh --refresh-data   first regenerate data/certificate_cases.json
#                               from an exact verifier run (about 92 s CPU)
#
# Needs uv (Python runs in a throwaway environment with matplotlib and numpy)
# and a TeX installation with latexmk (the figures use usetex).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
PY=(uv run --no-project --with matplotlib --with numpy python -I)

if [[ "${1:-}" == "--refresh-data" ]]; then
  out="$(mktemp)"
  uv run --no-project python -I "$REPO/verifier/verify_constant_gain_1_40.py" \
    "$REPO/certificate/constant-gain-1-40-K160-c60.json.gz" | tee "$out"
  "${PY[@]}" "$HERE/src/extract_data.py" certificate "$out"
  rm -f "$out"
fi

for name in fig_examples fig_circles_around_p fig_isosceles fig_annulus \
            fig_distance_lists fig_chord fig_bisector_gap fig_profile \
            fig_farkas fig_certificate; do
  echo "== $name"
  "${PY[@]}" "$HERE/src/$name.py"
done

echo "== fig_roadmap"
tmp="$(mktemp -d)"
if ! latexmk -pdf -cd -interaction=nonstopmode -halt-on-error \
     -outdir="$tmp" "$HERE/src/fig_roadmap.tex" > "$tmp/build.log" 2>&1; then
  cat "$tmp/build.log"
  exit 1
fi
cp "$tmp/fig_roadmap.pdf" "$HERE/fig_roadmap.pdf"
rm -rf "$tmp"
echo "done: $(ls "$HERE"/fig_*.pdf | wc -l | tr -d ' ') PDFs in $HERE"
