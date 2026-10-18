#!/bin/bash
# Reclaim space and inodes in a CSR run directory.
#
#   bash cleanup.sh /scratch/user/uquadhik/CSR_crystal_v2            # list only (default)
#   bash cleanup.sh /scratch/user/uquadhik/CSR_crystal_v2 --delete   # actually remove
#
# Lists by default and NEVER deletes without --delete. Read the listing before
# passing --delete.
#
# Both limits matter on Bunya: the GB quota and the file-count quota. ASH leaves
# dozens of small working files per optimisation, so a few hundred runs can eat
# the inode budget while using little space.
#
# What this NEVER touches:
#   preprocess_results/  md_results/  qmmm_results/  config.yaml
#   QM_MM.dcd  QM_MM_lastframe.pdb  QM_MM_optimized.pdb
#   scan_energies.csv  scan_*.pdb  nac_frame*.pdb  *.slurm  *.out  *.err
# i.e. every input, every converged structure and every result table survives.
set -uo pipefail

RUN=${1:-}
MODE=${2:-list}
if [ -z "$RUN" ] || [ ! -d "$RUN" ]; then
    echo "usage: bash cleanup.sh <run directory> [--delete]" >&2
    exit 1
fi
cd "$RUN" || exit 1

if [ "$MODE" = "--delete" ]; then
    echo "*** DELETING ***"
else
    echo "*** listing only; pass --delete to remove ***"
fi
echo

total_kb=0
total_n=0

# report <label> <find-expression...>
report() {
    local label=$1; shift
    local list kb n
    list=$(find "$@" 2>/dev/null) || true
    if [ -z "$list" ]; then
        printf '%-46s %10s %8s\n' "$label" "-" "0"
        return
    fi
    n=$(printf '%s\n' "$list" | wc -l)
    kb=$(printf '%s\n' "$list" | xargs -d '\n' -r du -sk 2>/dev/null | awk '{s+=$1} END{print s+0}')
    printf '%-46s %9s M %8s\n' "$label" "$((kb / 1024))" "$n"
    total_kb=$((total_kb + kb))
    total_n=$((total_n + n))
    if [ "$MODE" = "--delete" ]; then
        printf '%s\n' "$list" | xargs -d '\n' -r rm -rf --
    fi
}

printf '%-46s %11s %8s\n' "category" "size" "files"
printf '%-46s %11s %8s\n' "--------" "----" "-----"

# 1. ASH/geomeTRIC working files, regenerable, and the bulk of the file count.
#    Kept out of the protected directories by -path pruning.
report "ASH working files (xtb_, gradient, pcgrad...)" . \
    -path ./preprocess_results -prune -o -path ./md_results -prune -o \
    -path ./qmmm_results -prune -o \
    \( -name 'xtb_.*' -o -name 'pcgrad' -o -name 'gradient' -o -name 'energy' \
       -o -name 'charges' -o -name 'Fragment-currentgeo.*' \
       -o -name 'initialxyzfiletric.xyz' -o -name 'constraints.txt' \
       -o -name '*.engrad' -o -name 'ActiveRegion.xyz' \) -type f -print

# 2. geomeTRIC trajectories: large, and superseded by the final PDB once converged.
report "geomeTRIC optimisation trajectories" . \
    -path ./preprocess_results -prune -o -path ./md_results -prune -o \
    -path ./qmmm_results -prune -o \
    \( -name 'geometric_OPTtraj*.xyz' -o -name 'geometric_OPTtraj.tmp' \
       -o -name 'HT.v_*.txt' \) -type f -print

# 3. Runs that were killed or cancelled and produced nothing usable.
# Before removing the killed scans, keep their energy tables: a few kB each, and
# the only record of how far those runs got.
if [ -d scans ]; then
    mkdir -p scans/killed_runs_csv
    for d in scans/scan_H[0-9]_f*; do
        case "$d" in */scan_H1_f00497|*/scan_H2_f00395) continue ;; esac
        [ -f "$d/scan_energies.csv" ] || continue
        cp -n "$d/scan_energies.csv" "scans/killed_runs_csv/$(basename "$d").csv" 2>/dev/null
    done
fi
report "out-of-memory / cancelled 1D scans (CSVs kept)" \
    scans -maxdepth 1 -type d -name 'scan_H[0-9]_f*' \
    ! -name 'scan_H1_f00497' ! -name 'scan_H2_f00395' -print
report "cancelled 2D grids" scans -maxdepth 1 -type d -name 'scan2d_*' -print
report "abandoned GPU benchmark" protdft -maxdepth 1 -type d -name 'gpubench' -print
report "superseded QM/MM probe directories" . -maxdepth 1 -type d \
    \( -name 'qmmm_results_*cancelled' -o -name 'qmmm_results_*unbound' \) -print

# 4. Config backups: keep the three most recent.
old_baks=$(ls -1t config.yaml.bak* 2>/dev/null | tail -n +4)
if [ -n "$old_baks" ]; then
    kb=$(printf '%s\n' "$old_baks" | xargs -d '\n' -r du -sk | awk '{s+=$1} END{print s+0}')
    n=$(printf '%s\n' "$old_baks" | wc -l)
    printf '%-46s %9s M %8s\n' "config backups beyond the newest 3" "$((kb / 1024))" "$n"
    total_kb=$((total_kb + kb)); total_n=$((total_n + n))
    [ "$MODE" = "--delete" ] && printf '%s\n' "$old_baks" | xargs -d '\n' -r rm -f --
fi

echo
printf 'total %d MB in %d files\n' "$((total_kb / 1024))" "$total_n"
echo
echo "Not covered here, decide separately:"
echo "  - scans/scan_H1_f00497 and scans/scan_H2_f00395 keep their results;"
echo "    their working files are cleaned by categories 1-2 above."
echo "  - scans/frames_H*/nac_frame*.pdb are start structures: small, kept."
echo "  - qmmm_results/QM_MM.dcd is the 50 ps QM/MM MD trajectory: keep it,"
echo "    every NAC frame and future reanalysis comes from it."
echo "  - GPU python packages in the ASH conda env (cupy, nvidia-*-cu12,"
echo "    gpu4pyscf) are ~2-3 GB in HOME, not scratch. If GPU work is dropped:"
echo "      pip uninstall -y gpu4pyscf-cuda12x gpu4pyscf-libxc-cuda12x \\"
echo "          cupy-cuda12x fastrlock pyscf-dispersion \\"
echo "          \$(pip list --format=freeze | grep -o '^nvidia-[^=]*')"
echo "    Check numpy/openmm/mdtraj still import afterwards."
