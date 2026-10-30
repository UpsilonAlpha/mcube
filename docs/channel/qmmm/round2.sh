#!/bin/bash
# Second round of the QM/MM exchange energies (README "Second round"; decision log C22), one step
# per call so each can be checked before the next:
#
#   bash round2.sh seeds        Ca swapped into a bound La and a bound Lu frame (run dirs CA_la, CA_lu)
#   bash round2.sh mm           submit the two Ca windows (GPU) and the Ca bulk-water frames (CPU)
#   bash round2.sh lachbe       frames + submit: a third site source, La under the softer C4
#   bash round2.sh variants     D3 on every finished job, wB97X-V on two frames per source (copies)
#   bash round2.sh potential    far-field MM potential at the ion (writes qmmm/potential.json)
#   bash round2.sh complex      ion + capped Glu boxes for the calibration (CPU MM, 2 ns each)
#   bash round2.sh ca           frames + submit for the Ca sources, once `mm` has finished
#   bash round2.sh calib        frames + submit for the calibration boxes, once `complex` has finished
#   bash round2.sh inbox        the water reference inside the membrane box (PMF outer windows), Ca and La
set -e
S=/scratch/user/uquadhik
Q=$S/mcube/docs/channel/qmmm
SY=$S/mcube/docs/channel/system
B=$S/mcube/docs/channel/bulk
OUT=$S/channel/qmmm
CAL=$S/channel/qmmm_cal
IN=$S/channel/qmmm_inbox
BASIS="--basis-svp $S/channel/dft/def2-svp_lcecp1.nw --basis-tzvp $S/channel/dft/def2-tzvp_lcecp1.nw"
module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
conda activate ASH
export PYTHONPATH=$S/mcube/src PYTHONUNBUFFERED=1
P=$CONDA_PREFIX/bin/python
RES=(k=2.5 time=10000 'frames_at=[0.0]' traj_frequency=6250)

# submit the jobs of LIST from line FIRST on (lines already submitted stay out)
submit_from() {   # LIST FIRST NAME TIME
    local n=$(grep -c . "$1")
    [ "$n" -ge "$2" ] || { echo "no new jobs in $1"; return; }
    (cd $(dirname "$1") && sbatch -J $3 --array=$2-$n%10 -t $4 --export=ALL,LIST=$1 $Q/exchange.slurm)
}
lines() { [ -f "$1" ] && grep -c . "$1" || echo 0; }

case $1 in
seeds)
    $P $SY/swap_seed.py --window $S/channel/system/LA/umbrella_res_035 --frame 159 --ion CA \
        --base $S/channel/system/CA/config.yaml --out $S/channel/system/CA_la --set "${RES[@]}"
    $P $SY/swap_seed.py --window $S/channel/system/LU/umbrella_res_chbe_035 --frame 159 --ion CA \
        --base $S/channel/system/CA/config.yaml --out $S/channel/system/CA_lu --set "${RES[@]}" ;;
mm)
    for D in CA_la CA_lu; do
        (cd $S/channel/system/$D && sbatch -J res_$D --export=ALL,FIRST=35,LAST=35,SECTION=umbrella_res $SY/umbrella_windows.slurm)
    done
    (cd $S/channel/bulk && sbatch -J bulkf_CA --export=ALL,RUN=$S/channel/bulk/CA_1264_frames,NS=2,SKIP_QM=1,MMOPTS="--temp 310 --dcd-every 50" $B/bulk.slurm CA 1264) ;;
lachbe)
    first=$(( $(lines $OUT/jobs_site.txt) + 1 ))
    $P $Q/exchange.py frames --out $OUT --site LA@LAchbe=$S/channel/system/LA/umbrella_res_chbe_035 --actradius 0.5 $BASIS
    submit_from $OUT/jobs_site.txt $first qx_site 24:00:00 ;;
variants)
    $P $Q/exchange.py variant --out $OUT --dest $S/channel/qmmm_d3 --set sp_dispersion=d3bj
    $P $Q/exchange.py variant --out $OUT --dest $S/channel/qmmm_wb97xv --set sp_functional=wb97x-v --frames-per-source 2
    for V in qmmm_d3 qmmm_wb97xv; do
        submit_from $S/channel/$V/jobs_site.txt 1 v_${V#qmmm_}_s 8:00:00
        submit_from $S/channel/$V/jobs_water.txt 1 v_${V#qmmm_}_w 4:00:00
    done ;;
potential)
    $P $Q/exchange.py potential --out $OUT --pme ;;
complex)
    # the Glu closest to the ion in the first La site frame, in that frame's own numbering
    F=$OUT/site/LA_f0159/LA/frame.pdb
    GLU=$($P -c "
from openmm.app import PDBFile; import numpy as np
p = PDBFile('$F'); x = np.array(p.positions._value)
ion = next(a for a in p.topology.atoms() if a.residue.name == 'LA')
cd = [a for a in p.topology.atoms() if a.residue.name == 'GLU' and a.name == 'CD']
a = min(cd, key=lambda a: np.linalg.norm(x[a.index] - x[ion.index]))
print(f'{a.residue.chain.id}:{a.residue.id}')")
    echo "Glu for the complex: $GLU"
    for ION in CA LA LU; do
        (cd $S/channel/bulk && sbatch -J cplx_$ION --export=ALL,RUN=$S/channel/bulk/${ION}_glu_frames,NS=2,SKIP_QM=1,MMOPTS="--temp 310 --dcd-every 50 --ligand-from $F --ligand-resid $GLU" $B/bulk.slurm $ION 1264)
    done ;;
ca)
    first=$(( $(lines $OUT/jobs_site.txt) + 1 )); wfirst=$(( $(lines $OUT/jobs_water.txt) + 1 ))
    $P $Q/exchange.py frames --out $OUT --site CA@CAla=$S/channel/system/CA_la/umbrella_res_035 \
        --site CA@CAlu=$S/channel/system/CA_lu/umbrella_res_035 --water CA=$S/channel/bulk/CA_1264_frames \
        --actradius 0.5 $BASIS
    submit_from $OUT/jobs_site.txt $first qx_site 24:00:00
    submit_from $OUT/jobs_water.txt $wfirst qx_water 8:00:00 ;;
calib)
    $P $Q/exchange.py frames --out $CAL --complex CA@acCA=$S/channel/bulk/CA_glu_frames \
        --complex LA@acLA=$S/channel/bulk/LA_glu_frames --complex LU@acLU=$S/channel/bulk/LU_glu_frames \
        --actradius 0.5 $BASIS
    # the water reference is the one the site uses (its jobs stay in $OUT): link its entries in
    $P -c "
import json; a = json.load(open('$CAL/frames.json')); b = json.load(open('$OUT/frames.json'))
a.update({k: v for k, v in b.items() if v['kind'] == 'water'}); json.dump(a, open('$CAL/frames.json', 'w'), indent=1)"
    submit_from $CAL/jobs_site.txt 1 qc_site 8:00:00 ;;
inbox)
    # The water reference inside the membrane box (C22): the last frames of the La and Ca PMFs'
    # outer windows, ion in the bulk water layer, Ca and La swapped into each. The site jobs stay
    # in $OUT; their entries are linked in (rerun after `ca` to take in the Ca site sources).
    first=$(( $(lines $IN/jobs_water.txt) + 1 ))
    U=$S/channel/system
    $P $Q/exchange.py frames --out $IN --ions CA LA --actradius 0.5 $BASIS \
        --bulk LA@LAbulk=$U/LA/umbrella_010:$U/LA/umbrella_011:$U/LA/umbrella_062:$U/LA/umbrella_063:$U/LA/umbrella_064 \
        --bulk CA@CAbulk=$U/CA/umbrella_000:$U/CA/umbrella_001:$U/CA/umbrella_002:$U/CA/umbrella_069:$U/CA/umbrella_070
    $P -c "
import json; a = json.load(open('$IN/frames.json')); b = json.load(open('$OUT/frames.json'))
a.update({k: v for k, v in b.items() if v['kind'] == 'site'}); json.dump(a, open('$IN/frames.json', 'w'), indent=1)"
    submit_from $IN/jobs_water.txt $first qx_bulk 24:00:00 ;;
*)
    sed -n 2,15p $0; exit 1 ;;
esac
