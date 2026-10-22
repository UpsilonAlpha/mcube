#!/bin/bash
#SBATCH --partition=gpu_rocm --qos=gpu --gres=gpu:mi355x:1 --account=a_imb_ccbcd
#SBATCH --cpus-per-task=8 --mem=96G --time=12:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# AlphaFold 2.3.2 through RCC's `alphafold_custom` (nothing preselected), with every
# database given explicitly. Same paths as the module's own `alphafold` /
# `alphafold_multimer` functions, except pdb_seqres: both modules point
# `alphafold_multimer` at .../databses/pdb_seqres/... (missing "a"), which does not exist.
#
# Default: an MI355X on gpu_rocm, with the alphafold/2.3.2-rocm build. To run on NVIDIA
# instead, override the partition and GPU on the sbatch line; the script picks the CUDA
# build from the partition. Only A100 there: RCC's CUDA container wrapper notes it is
# incompatible with H100 (CUDA < 11.8).
#
#   sbatch -J af2_mono   docs/channel/af2/af2.slurm monomer  docs/channel/af2/mature.fasta        $OUT/monomer
#   sbatch -J af2_trimer docs/channel/af2/af2.slurm multimer docs/channel/af2/mature_trimer.fasta $OUT/trimer
#   sbatch -p gpu_cuda --gres=gpu:a100:1 -J ... (same arguments)          # the CUDA build
#
# Monomer uses monomer_ptm (pTM and PAE as well as pLDDT). Multimer makes one prediction
# per model (5 models); raise NPRED for the full 25.

PRESET=$1
FASTA=$(readlink -f "$2")
OUT=$3
NPRED=${NPRED:-1}
mkdir -p "$OUT"

# The modules define alphafold_custom as a shell function; each branch defines the same
# thing if a non-interactive batch shell did not receive it
if [ "$SLURM_JOB_PARTITION" = "gpu_rocm" ]; then
    module load alphafold/2.3.2-rocm
    export APPTAINERENV_PREPEND_PATH=${APPTAINERENV_PREPEND_PATH:-/opt/hh-suite/bin}
    type alphafold_custom >/dev/null 2>&1 || alphafold_custom() {
        apptainer exec --rocm /sw/containers/local/rcc/alphafold2-amd-gpu_rocm6.1.1.sif \
            python /opt/alphafold/run_alphafold.py "$@"; }
    GPUINFO=$(rocm-smi --showproductname 2>/dev/null | grep -m1 -i 'card series' || echo rocm)
else
    module load alphafold/2.3.2
    type alphafold_custom >/dev/null 2>&1 || alphafold_custom() {
        /sw/containers/local/rcc/CVL/bin/cvl-alphafold_2.3.2-noninteractive.sh "$@"; }
    GPUINFO=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null)
fi
D=/scratch/opendata/protein/AlphaFold/databases

COMMON="--data_dir=$D --use_gpu_relax=true --max_template_date=2024-05-05
 --uniref90_database_path=$D/uniref90/uniref90.fasta
 --mgnify_database_path=$D/mgnify/mgy_clusters_2018_12.fa
 --template_mmcif_dir=$D/pdb_mmcif/mmcif_files/
 --obsolete_pdbs_path=$D/pdb_mmcif/obsolete.dat
 --bfd_database_path=$D/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
 --uniref30_database_path=$D/uniclust30/uniclust30_2018_08/uniclust30_2018_08
 --fasta_paths=$FASTA --output_dir=$OUT"

case "$PRESET" in
    monomer)
        EXTRA="--model_preset=monomer_ptm --pdb70_database_path=$D/pdb70/pdb70" ;;
    multimer)
        EXTRA="--model_preset=multimer --num_multimer_predictions_per_model=$NPRED
 --uniprot_database_path=$D/uniprot/uniprot.fasta
 --pdb_seqres_database_path=$D/pdb_seqres/pdb_seqres.txt" ;;
    *)
        echo "usage: af2.slurm monomer|multimer FASTA OUTDIR" >&2; exit 2 ;;
esac

echo "host $(hostname)  partition $SLURM_JOB_PARTITION  gpu $GPUINFO"
echo "preset $PRESET  fasta $FASTA  out $OUT"
alphafold_custom $COMMON $EXTRA
RC=$?
echo "RC=$RC"
exit $RC
