# CSR-1 / meropenem: cluster scripts

Four files, not one per calculation. Everything specific to this system lives in
`site.yaml`; the scripts stay general, so a new calculation is a YAML edit and a
new system is a new YAML file.

| File | What it is |
| --- | --- |
| `site.yaml` | the active site (atoms, candidate bases, near-attack criteria) and every calculation we run |
| `sections.py` | adds calculations from `site.yaml` to `config.yaml`, refusing to write if a selector does not resolve |
| `site.py` | all the geometry analysis: `sphere`, `nucleophiles`, `acceptors`, `joint` |
| `run.slurm` | one job script for every mcube stage; resources come from the `sbatch` line |
| `cleanup.sh` | reclaim space and inodes in a run directory |

They are kept in the repo and pulled on Bunya rather than typed into the browser
terminal, which silently truncates long payloads.

## Running something

```bash
cd /scratch/user/uquadhik/mcube && git pull -q
cd /scratch/user/uquadhik/CSR_crystal_v2
cp config.yaml config.yaml.bak          # the stages write results back into it

D=/scratch/user/uquadhik/mcube/docs/csr
python $D/sections.py $D/site.yaml --list          # what is defined
python $D/sections.py $D/site.yaml scan_znco       # add it, checking selectors
sbatch -J znco --cpus-per-task=8 --mem=16G $D/run.slurm \
       xtb scan scan_znco protdft/zncheck/relax/QM_MM_optimized.pdb
```

`sections.py` writes `<stage>_<name>`, so the `scan_znco` entry in `site.yaml`
becomes section `scan_znco`, and `cleanR_lysH` becomes `qmmm_cleanR_lysH`.

## Choosing resources

The full measured table is in
[docs/running-on-bunya.md](../running-on-bunya.md#resources-by-job-type); the
short version, which is why `run.slurm` takes a profile:

| Work | CPUs | CPU efficiency | Memory used |
| --- | --- | --- | --- |
| DFT single points | 16 | **96%** | 56 GB of 96 |
| xTB optimisation or scan | 8 | **94%** (geomeTRIC, not xTB) | ~1 GB |
| xTB/MM MD | 8 | **30%** | 0.5 GB of 32 |

DFT scales, so give it cores — 64 is available and worth using, and it is the
one stage that needs real memory. xTB does not scale: it saturates around four
threads, so make xTB work finish sooner by running points in parallel rather
than adding threads.

The `xtb` profile binds threads (worth 2–9× on SMT nodes, where allocated CPUs
are half physical cores). The `dft` profile deliberately does not, and sets
`OMP_NUM_THREADS=1`, because PySCF runs in-process, sets its own thread count
from `dft_cores`/`sp_cores`, and would otherwise end up pinned to one core.

## Analysis

```bash
D=/scratch/user/uquadhik/mcube/docs/csr
python $D/site.py sphere       --site=$D/site.yaml structure.pdb
python $D/site.py nucleophiles --site=$D/site.yaml structure.pdb
python $D/site.py acceptors    --site=$D/site.yaml traj.dcd topology.pdb
python $D/site.py joint        --site=$D/site.yaml traj.dcd topology.pdb
```

- **`sphere`** — what is coordinated to the metal, the monitored distances, and
  how many hydrogens sit on each candidate oxygen and on the lysine. The quickest
  way to see whether a structure is what you think it is.
- **`nucleophiles`** — for one structure: every water near the electrophile, its
  attack distance and angle, whether it is metal-bound, and for each of its
  protons what lies within hydrogen-bond reach.
- **`acceptors`** — over a trajectory, during the frames where a water is
  genuinely on the Bürgi–Dunitz line, what is near its protons *at all*. Use it
  when you do not want to presume which base is involved.
- **`joint`** — over a trajectory, does any single water both reach the attack
  geometry *and* donate to a named base, directly or through one bridging water?
  A general-base mechanism needs both at once, and this is the test that ruled
  the direct Asp92 route out twice.

## Why the selector check exists

`sections.py` resolves every `chain:resid:atom` in a section against that
section's own start structure and refuses to write if any is missing. Two days
of compute were lost to a 2D grid whose second coordinate began 4.35 Å apart —
which drags a group across space rather than transferring a proton. Check the
coordinate against the structure it will actually start from.

## Email notifications

`run.slurm` carries `#SBATCH --mail-type=END,FAIL` but deliberately **no
`--mail-user`**: this repository is public, so an address committed here would
be scraped. Pass one at submit time, which overrides the directives:

```bash
sbatch --mail-user=you@example.com ... run.slurm ...
```

For a job already queued or running: `scontrol update JobId=<id>
MailUser=you@example.com MailType=END,FAIL`.
