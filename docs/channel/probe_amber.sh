#!/bin/bash
# Read-only: find AmberTools on Bunya and print the 12-6-4 rows mcube needs, so the values
# vendored in src/mcube/ions.py can be checked against them.
#
#   bash docs/channel/probe_amber.sh > $L/probe_amber.log 2>&1
#
# Prints: the Amber modules; for each install, the La/Lu/Y/Gd rows of
# frcmod.ions34lsm_1264_tip3p, the Ca/Mg/Zn rows of frcmod.ions234lm_1264_tip3p (the +2 set
# mcube lacks), and the sha256 of lj_1264_pol.dat.

source /etc/profile.d/modules.sh 2>/dev/null || true
echo "== modules"
module -t avail 2>&1 | grep -i -E '^amber' || echo "no amber module"

roots=""
for m in $(module -t avail 2>&1 | grep -i -E '^amber'); do
    p=$(module show "$m" 2>&1 | grep -E 'AMBERHOME' | grep -oE '/[^ ")]+' | head -1)
    [ -n "$p" ] && roots="$roots $p"
done
[ -n "$AMBERHOME" ] && roots="$roots $AMBERHOME"
echo "== roots:$roots"

for r in $roots; do
    for f in $(find "$r" -maxdepth 6 \( -name 'frcmod.ions34lsm_1264_tip3p' \
            -o -name 'frcmod.ions234lm_1264_tip3p' -o -name 'lj_1264_pol.dat' \) 2>/dev/null); do
        echo "== $f"
        case "$f" in
            *34lsm_1264_tip3p) grep -E '^(La3\+|Lu3\+|Y3\+|Gd3\+) ' "$f" ;;
            *234lm_1264_tip3p) head -1 "$f"; grep -E '^(Ca2\+|Mg2\+|Zn2\+) ' "$f" ;;
            *lj_1264_pol.dat) sha256sum "$f"; grep -E '^(OW|CT|Ca2\+) ' "$f" ;;
        esac
    done
done
echo "== done"
