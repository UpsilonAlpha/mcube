#!/usr/bin/env bash
# Pull the CSR coordinate files worth having locally.
#
# RUN THIS YOURSELF. There is no SSH private key on this machine - ~/.ssh/ holds
# only known_hosts - so scp has to ask for your Bunya credentials, and that is
# the one step that cannot be driven from the assistant's side. Everything below
# is a plain scp; nothing is written on Bunya.
#
#   bash docs/csr/structures/pull-from-bunya.sh
#
# Files land next to this script. *.pdb here is gitignored on purpose.
set -euo pipefail

HOST=${HOST:-uquadhik@bunya.rcc.uq.edu.au}
S=/scratch/user/uquadhik
DEST=$(cd "$(dirname "$0")" && pwd)

echo "pulling into $DEST from $HOST"

# One scp invocation per group, so one authentication covers the group and a
# missing file does not abort the rest.

# --- the product arm, restarted 2026-11-09 (D249) ---
scp "$HOST:$S/CSR/CSR_prod/md_prodpull_results/NVTtrajectory_lastframe.pdb" \
    "$DEST/prodpull_lastframe.pdb"
scp "$HOST:$S/CSR/CSR_prod/preprocess_results/finalsystem.pdb" \
    "$DEST/prod_finalsystem.pdb"

# --- the NaN evidence (D249): stable to the last frame, then died ---
scp "$HOST:$S/CSR/CSR_prod/md_prod100c_results/nvtsim.csv" \
    "$DEST/prod100c_nvtsim.csv"
scp "$HOST:$S/CSR/CSR_prod/md_prodpull_results/nvtsim.csv" \
    "$DEST/prodpull_nvtsim.csv"

# --- the reactant arm: the frame every QM/MM result descends from ---
scp "$HOST:$S/CSR/CSR_mich/md_relax10_results/NVTtrajectory_lastframe.pdb" \
    "$DEST/relax10_lastframe.pdb"

# --- the two converged QM/MM minima (D240, D244) ---
scp "$HOST:$S/CSR/CSR_mich/protdft/qmmm_aspxtb/QM_MM_optimized.pdb" \
    "$DEST/aspxtb_minimum.pdb"
scp "$HOST:$S/CSR/CSR_mich/protdft/qmmm_aspxtb8/QM_MM_optimized.pdb" \
    "$DEST/aspxtb8_minimum.pdb"

# --- the def2-SVP optimisation trajectory, ~42 MB (D237/D243) ---
# The wall writes no QM_MM_optimized.pdb; this is what traj_frame.py recovers
# the final geometry from. Comment out if you only want the recovered frame.
scp "$HOST:$S/CSR/CSR_mich/protdft/qmmm_aspdft/geometric_OPTtraj_Full.xyz" \
    "$DEST/aspdft_OPTtraj_Full.xyz"

echo "done. both arms live under \$S/CSR/ - CSR_prod and CSR_mich."
