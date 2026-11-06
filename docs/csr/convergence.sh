#!/bin/bash
# Run convergence.py on every trajectory that bears on the 10 ns vs 100 ns
# question (D225). Read-only: writes nothing but its own stdout.
#
#   wt100c    100 ns  - condemned for SITE geometry (D194), but the global
#                       backbone convergence it reports is still valid, which is
#                       the whole reason we need not run a fresh 100 ns to ask
#   wt100b    100 ns  - the other long one (condemned, D177); a second opinion
#   eq6        10 ns  - the equilibration this project actually did
#   asppull3    2 ns  - the repaired structure; has ITS site finished moving?
set -uo pipefail

R=/scratch/user/uquadhik/mcube/docs/csr
S=/scratch/user/uquadhik/CSR

module load miniforge/26.1.0-0
# shellcheck disable=SC1091
source "$ROOTMINIFORGE/etc/profile.d/conda.sh"
conda activate ASH

echo "=== date: $(date)"
cd "$S"
for d in CSR_wt/md_wt100c_results \
         CSR_wt/md_wt100b_results \
         CSR_mich/md_eq6_results \
         CSR_mich/md_asppull3_results; do
  if [ -d "$d" ]; then
    python "$R/convergence.py" "$d"
  else
    echo "=== $d: absent"
  fi
done
echo "=== ALLDONE"
