# CSR-1 / meropenem: cluster scripts

Everything specific to this system lives in `site.yaml`; the scripts stay
general, so a new calculation is a YAML edit and a new system is a new YAML
file.

**The core four:**

| File | What it is |
| --- | --- |
| `site.yaml` | the active site (atoms, candidate bases, near-attack criteria) and every calculation we run |
| `sections.py` | adds calculations from `site.yaml` to `config.yaml`, refusing to write if a selector does not resolve |
| `site.py` | all the geometry analysis: `sphere`, `nucleophiles`, `acceptors`, `joint`, `poses` |
| `run.slurm` | one job script for every mcube stage; resources come from the `sbatch` line |
| `cleanup.sh` | reclaim space and inodes in a run directory |

**Measuring what a run produced:**

| File | What it is |
| --- | --- |
| `traj_frame.py` | writes one frame of an xyz trajectory onto a PDB template, validating the atom count and the element-vs-name order |
| `traj_contacts.py` | named active-site contacts over a trajectory - min, max, mean and the fraction of frames inside a cutoff, with the closest atom pair named (D111) |

**Building the ohH arm of the protonation test** (D112-D114). Only needed while
that comparison is open; see `site.yaml`'s `ohHmin` and `ohH` sections.

| File | What it is |
| --- | --- |
| `make_water_nuc.py` | turns the hydroxide into a water and Lys226 neutral, in an equilibrated frame: adds H04, renames OH to NUW, removes HZ3 with `--lys A:226`. `--floor` sets the clearance the new hydrogens need |
| `water_nuc.xml` | the water template. Declares nothing but the residue, borrowing TIP3P's own atom types, so it is parameter-identical to HOH - a template with its own types matches every bulk water by bond graph and OpenMM then refuses the system |
| `check_ohH_frame.py` | builds **both** arms' systems in full and asserts they differ by one hydrogen: same atom count, same total charge, nucleophile -1 to 0, lysine +1 to 0 |
| `check_ohH_system.py` | the same checks for a structure that has not been solvated yet, read from the parsed templates rather than a built system |
| `make_ohH_config.py` | derives an ohH config from the lysH one, for the abandoned fresh-preprocess route (D113). Kept because it is the only record of which three preprocess keys differ |

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
sbatch -J znco --cpus-per-task=8 --mem=16G $D/run.slurm scan scan_znco
```

The starting structure belongs in the section's `pdb_file` (a NEB's two
endpoints in `reactant_pdb` and `product_pdb`), and there is no longer any way to
pass one on the command line: `run.slurm` takes no third argument and no stage
takes `--pdb`. `md` had always read its structure from the config and ignored a
`--pdb`, which could silently equilibrate preprocess's own output instead (D114);
every stage now works that way, so `config.yaml` describes the run on its own.

Paths in `inputs` are read against the config's directory, whichever directory
the job script runs the stage from — so `md_results/lastframe.pdb` means the same
file to `md` and to `qmmm`. (`forcefields` is the exception: OpenMM resolves
those itself, from the directory the stage runs in, which is why they carry
`../`.) The active region is `actatoms` (index ranges, e.g. `4203-4634,5001`) or
`actatoms_from: <section>` to reuse the region another section recorded, in place
of the `active_atoms` file that used to be passed with `--actatoms`. A recorded
region carries the structure it was built on and a per-atom fingerprint, and
`actatoms_from` refuses a structure it does not describe — the D115 failure mode,
where indices saved against one atom order were applied to another.


`sections.py` writes `<stage>_<name>`, so the `scan_znco` entry in `site.yaml`
becomes section `scan_znco`, and `cleanR_lysH` becomes `qmmm_cleanR_lysH`.

## Choosing resources

Only `--cpus-per-task`, `--mem` and `--time` are yours to pick. Thread placement
is derived from the stage and the config and applied by `run.slurm`, which also
warns when the allocation is bigger than the stage can use — so there is nothing
else to get right. `run.slurm`'s header carries a worked `sbatch` line per stage.

| Stage | CPUs | Memory |
| --- | --- | --- |
| `preprocess` | 2 | 4G |
| `md` | 16 | 16G |
| `qmmm`, `scan`, `ts` | 8 | 16G |
| `neb` | 16 | 32G |
| `sp` | 32 | 96G |

The measurements, and what dominates each stage, are in
[docs/running-on-bunya.md](../running-on-bunya.md#resources-by-job-type).

## Analysis

```bash
D=/scratch/user/uquadhik/mcube/docs/csr
python $D/site.py sphere       --site=$D/site.yaml structure.pdb
python $D/site.py nucleophiles --site=$D/site.yaml structure.pdb
python $D/site.py acceptors    --site=$D/site.yaml traj.dcd topology.pdb
python $D/site.py joint        --site=$D/site.yaml traj.dcd topology.pdb
python $D/site.py poses        --site=$D/site.yaml docked.pdbqt receptor.pdb [ligand.pdb]
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
- **`poses`** — the docked poses in a PDBQT, each one's Vina score and how far
  its centroid sits from the metal, and — given the ligand structure actually in
  use — which pose that structure is. Vina scores within a couple of kcal/mol of
  each other do not pick a pose, so check what the alternatives look like before
  building a mechanism on one of them.

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
