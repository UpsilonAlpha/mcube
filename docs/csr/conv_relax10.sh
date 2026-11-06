#!/bin/bash
# Did `relax10` converge? Read-only: writes nothing but its own stdout.
#
# The gate before QM/MM (D228): require eq6's signature - CA RMSD flat over the
# last five deciles and NVT box-volume drift under ~1 nm3. `asppull3`'s 2 ns
# drifted +5.396 nm3 against eq6's +0.653, which is why this run exists. A drift
# still near +5.4 after 10 ns is the restraint set fighting the solvent, not a
# run needing longer, and in that case QM/MM does NOT start.
#
# eq6 is re-run alongside it so the comparison is against numbers produced by
# the same code on the same day rather than against a figure copied out of a log.
set -uo pipefail

R=/scratch/user/uquadhik/mcube/docs/csr
S=/scratch/user/uquadhik/CSR

module load miniforge/26.1.0-0
# shellcheck disable=SC1091
source "$ROOTMINIFORGE/etc/profile.d/conda.sh"
conda activate ASH

echo "=== date: $(date)"
cd "$S"
for d in CSR_mich/md_relax10_results CSR_mich/md_eq6_results; do
  if [ -d "$d" ]; then
    python "$R/convergence.py" "$d"
  else
    echo "=== $d: absent"
  fi
done
echo "=== ALLDONE"
