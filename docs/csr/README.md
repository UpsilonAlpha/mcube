# CSR-1 job scripts

Project-specific scripts for the CSR-1 / meropenem study (see `../csr-decision-log.md`).
They are kept here so the cluster can `git pull` them instead of having them retyped
into a browser terminal, which corrupts long payloads.

Run from the run directory, e.g. `/scratch/user/uquadhik/CSR_crystal_v2`:

```bash
cp config.yaml config.yaml.bak            # keep a backup
python /scratch/user/uquadhik/mcube/docs/csr/mk_cleanR.py
cp /scratch/user/uquadhik/mcube/docs/csr/CleanR.slurm .
sbatch CleanR.slurm
```

## `mk_cleanR.py` + `CleanR.slurm` — D42, the clean four-coordinate proton check

Adds the sections `qmmm_cleanR_lysH` and `qmmm_cleanR_ohH` and prints the resolved
atom index for every constrained selector (anything printed as `MISSING` means the
section is wrong — do not submit).

Both states start from `md_results/NVTtrajectory_lastframe.pdb`, the pre-QM/MM-MD
NVT end point, which already holds the canonical reactant: Zn + His88 + His164 +
OH- (O01, one H) + H2O (HOH485), with Lys226 as NH3+ and Asp92 OD2 at its
crystallographic 2.79 A.

The same two sphere constraints go on **both** states, so the energy difference
reflects only where the proton sits:

| Constraint | Value | Why |
| --- | --- | --- |
| Zn - Asp92 OD2 | 2.80 A | crystallographic; stops the migration to 2.07-2.11 A seen in the QM/MM MD |
| Zn - mer O2x | 2.90 A | carbonyl polarised, not inner-sphere; the preprocess metal restraint had pinned it at 1.94 A, shorter than the Zn-OH- bond |

Each state then pins the shared proton on one side, so neither can collapse into
the other:

| Section | Pinned | Gives |
| --- | --- | --- |
| `qmmm_cleanR_lysH` | NZ(Lys226)-HZ1 = 1.03 A | proton on the lysine; O01 stays hydroxide |
| `qmmm_cleanR_ohH` | O01-HZ1 = 0.98 A | proton on O01; O01 becomes water |

Pinning both sides is a change from D30, where `ohH` was left unconstrained and
relied on the transfer happening by itself. In a clean sphere it may not, which
would have made the two states identical and the comparison empty.

`CleanR.slurm` runs the two sequentially (16 CPUs, 96 GB, 2 days) and gives `ohH`
lysH's `active_atoms`, so both use the same active region. Read the result as
`sp_energy_disp` for each section: dE = E(ohH) - E(lysH), negative meaning the
transfer to the hydroxide is still favourable.

## `mk_gpubench.py` + `GpuBench.slurm` — D45, DFT/MM gradient timing on a GPU

Needs `gpu4pyscf` in the ASH env (installed 2026-10-18 with `cupy-cuda12x==13.3.0`
and `numpy==1.26.4` pinned, so the install was purely additive and numpy stayed at
1.26.4 — a newer cupy wants numpy >= 2, which would break OpenMM 8.4 and mdtraj).

Three levels, each a few constrained optimisation steps with no single point, so the
per-step wall time in `optimization_energies.log` is the DFT/MM gradient cost:
wB97X (comparable to the measured CPU cost), r2SCAN-D4 (the proposed geometry level),
and r2SCAN-D4 with `GPU_pcgrad` to price moving the MM point-charge gradient onto the
device. Requests one GPU on `gpu_cuda`, QOS `gpu`.

## Email notifications

Both job scripts carry `#SBATCH --mail-type=END,FAIL` but deliberately **no
`--mail-user`**: this repository is public, so an address committed here would be
scraped. SLURM mails the submitting user's configured address by default; if that
is not set up, pass one at submit time, which overrides the directives:

```bash
sbatch --mail-user=you@example.com CleanR.slurm
```

The working copies under `/scratch/.../CSR_crystal_v2/` are not in git, so an
address can be added there instead.

For a job that is already queued or running, change it in place:

```bash
scontrol update JobId=<id> MailUser=you@example.com MailType=END,FAIL
```
