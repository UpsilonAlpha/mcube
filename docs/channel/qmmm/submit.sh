#!/bin/bash
# Submit the QM/MM vertical-swap jobs `exchange.py frames` listed, as one array per kind.
#
#   bash submit.sh $S/channel/qmmm            # site and water
#   bash submit.sh $S/channel/qmmm site 10    # site only, at most 10 running at once
#
# Finished jobs exit at once, so resubmitting after a timeout only reruns what is missing.
OUT=${1:?usage: submit.sh OUT_DIR [site|water|both] [max running]}
KIND=${2:-both}
MAX=${3:-10}
HERE=$(cd "$(dirname "$0")" && pwd)
for K in site water; do
    [ "$KIND" = both ] || [ "$KIND" = "$K" ] || continue
    LIST=$OUT/jobs_$K.txt
    [ -s "$LIST" ] || { echo "no $LIST"; continue; }
    N=$(grep -c . "$LIST")
    # water boxes are ~4k atoms: much quicker than the ~64k-atom membrane frames
    TIME=24:00:00; [ $K = water ] && TIME=8:00:00
    (cd "$OUT" && sbatch -J qx_$K --array=1-$N%$MAX -t $TIME --export=ALL,LIST=$LIST "$HERE/exchange.slurm")
done
