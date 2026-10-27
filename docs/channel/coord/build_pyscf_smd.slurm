#!/bin/bash
#SBATCH --partition=general --account=a_imb_ccbcd --cpus-per-task=16 --mem=32G --time=3:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# PySCF with the SMD non-electrostatic (CDS) library, for `site_exchange.py --solvation smd`.
# The pip wheel is built without it ("compile with -DENABLE_SMD=ON"). This builds the same
# version as the ASH env from source into a venv on top of ASH (system site packages), so
# everything else (geomeTRIC, numpy) is shared:
#   sbatch -J build_smd docs/channel/coord/build_pyscf_smd.slurm
# then run with SMD=1 (site.slurm).

S=/scratch/user/uquadhik
V=$S/envs/smd
ASHPY=/home/uquadhik/.conda/envs/ASH/bin/python
VER=$($ASHPY -c 'import pyscf; print(pyscf.__version__)')
module load foss/2023a cmake/3.26.3-gcccore-12.3.0 || exit 1
[ -x $V/bin/python ] || $ASHPY -m venv --system-site-packages $V || exit 1
export CMAKE_CONFIGURE_ARGS="-DENABLE_SMD=ON" CMAKE_BUILD_PARALLEL_LEVEL=$SLURM_CPUS_PER_TASK
# --ignore-installed: otherwise pip sees ASH's wheel through the system site packages and skips
$V/bin/pip install -v --no-deps --ignore-installed --no-binary pyscf "pyscf==$VER" || exit 1
cd /tmp && $V/bin/python - <<'EOF'
import pyscf
from pyscf import gto, dft
from pyscf.solvent import smd
print("pyscf", pyscf.__version__, pyscf.__file__)
print("libsolvent", smd.libsolvent)
m = gto.M(atom="O 0 0 0; H 0 0 0.96; H 0.93 0 -0.24", basis="def2-svp", verbose=0)
for name in ("water", "diethylether", "propanonitrile"):
    mf = dft.RKS(m, xc="r2scan").density_fit().SMD()
    mf.with_solvent.solvent = name
    e = mf.kernel()
    print(name, smd.solvent_db[mf.with_solvent.solvent][:3], e)
print("La/Lu default SMD radii (A):", smd.smd_radii(0.0)[57], smd.smd_radii(0.0)[71])
EOF
echo "RC=$?"
