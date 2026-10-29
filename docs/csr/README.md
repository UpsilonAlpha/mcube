# CSR-1 / meropenem: cluster scripts

Everything specific to this system lives in `site.yaml`; the scripts stay
general, so a new calculation is a YAML edit and a new system is a new YAML
file.

**The core five:**

| File | What it is |
| --- | --- |
| `site.yaml` | the active site (atoms, candidate bases, near-attack criteria), the named Michaelis hypotheses, and every calculation we run |
| `sections.py` | adds calculations from `site.yaml` to `config.yaml`, refusing to write if a selector does not resolve — restraint pairs included (D128) |
| `site.py` | all the geometry analysis and the acceptance gate: `sphere`, `shell`, `nucleophiles`, `acceptors`, `joint`, `poses`, `graft`, `check` |
| `run.slurm` | one job script for every mcube stage; resources come from the `sbatch` line |
| `cleanup.sh` | reclaim space and inodes in a run directory |

**Reference values:**

| File | What it is |
| --- | --- |
| `refstruct.py` | `shell` and `centre` modes against a downloaded PDB, so every reference number in `site.yaml` is re-derivable rather than quoted (D125) |
| `reference-structures.md` | the PDB IDs, ligand codes, exact commands and measured tables behind those numbers |

## Read this before setting up a metal run

**The coordination detection cutoff is the most expensive parameter in this
project.** `preprocess.inputs.cutoff` decides which contacts become 50
kcal/mol/Å² restraints for the whole run, and a restraint only preserves a bond
if the contact *is* one. A Zn–O/N bond is 1.9–2.2 Å.

At `cutoff: 3.0` the crystal's second-shell Asp92 OD2 contact at **2.80 Å** was
detected as a ligand and held there; the carboxylate rotated its other oxygen
onto the metal at 1.99 Å, and two QM/MM arms ran on a five-coordinate zinc
against the crystal's four. **2.5 is defensible, 3.0 is not** (D120). The
library default is 2.8, so set it explicitly.

`md` reports every pair it is about to restrain that is longer than `LIGAND_MAX`
(2.5 Å), detected *and* configured, and records them in the section's
`coordination_beyond_bond_length`. Raise the threshold with
`preprocess.inputs.coordination_warn_above` where the chemistry warrants it — an
Ln–O bond reaches 2.4–2.5 Å. **Do not silence it by raising the detection
cutoff**, and do not hand-write a long restraint to correct a shell that is
already wrong: rebuild from a structure whose sphere is right (D130).

Holding a contact where the structure already has it is the rule this pipeline
runs on and it is sound. What it cannot do is decide whether the contact is a
bond — that is what the threshold above is for.

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

**One-shot, kept only as a record:** `migrate_paths.py` rewrote every config's
paths to be config-relative (D117) and has nothing left to migrate;
`make_ohH_config.py` describes an abandoned route. Neither is part of a current
workflow — do not reach for them when setting up a run.

## Running something

```bash
cd /scratch/user/uquadhik/mcube && git pull -q
cd /scratch/user/uquadhik/CSR_nac                # the live run directory
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
`../`.) The active region is not passed at all any more: it is derived from
`actradius` and `actregion_origin` on the structure the section runs on, and the
list it used is written to `active_atoms_used` beside the run, with its path,
count, structure and a per-atom fingerprint recorded in that section's results.
An `actatoms` or `actatoms_from` key is refused. That is the D115 lesson applied
to the one remaining index list: a region derived on one frame and applied to
another freezes the wrong atoms, and re-deriving is cheaper than guarding it.


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
python $D/site.py shell        --site=$D/site.yaml traj.dcd topology.pdb
python $D/site.py nucleophiles --site=$D/site.yaml structure.pdb
python $D/site.py acceptors    --site=$D/site.yaml traj.dcd topology.pdb
python $D/site.py joint        --site=$D/site.yaml traj.dcd topology.pdb
python $D/site.py poses        --site=$D/site.yaml docked.pdbqt receptor.pdb [ligand.pdb]
python $D/site.py graft        --site=$D/site.yaml docked.pdbqt ref.pdb receptor.pdb system.pdb <n> [--ref-frame=system]
python $D/site.py check        --site=$D/site.yaml structure.pdb --at=start|reactant --hypothesis=<name>
```

- **`sphere`** — what is coordinated to the metal, the monitored distances, and
  how many hydrogens sit on each candidate oxygen and on the lysine. The quickest
  way to see whether a structure is what you think it is. **One frame, so it
  cannot tell a ligand from a passer-by** — use `shell` on a trajectory for that.
- **`shell`** — the coordination sphere over a trajectory, ligand by ligand:
  each donor's **occupancy** (the fraction of frames it spends inside the
  cutoff) and the per-frame coordination number as a *distribution*, not an
  average, because a site flipping between 4 and 5 and a site sitting at 4.5 are
  different things. Donors just outside the cutoff are listed too, so a
  borderline call looks borderline instead of absent. This is the mode that
  found D120, and it is the right instrument for any claim about the sphere.
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
- **`graft`** — puts a chosen docked pose into an equilibrated system, mapping
  the pose onto the ligand structure already in use and dropping clashing
  waters. `--ref-frame=system` when the reference copy lives in the system
  rather than in the docking frame. Writes `<system>_pose<N>.pdb`; it does not
  equilibrate anything, and the result still has to pass `check`.
- **`check`** — the acceptance gate: accept or reject a structure against one
  **named** Michaelis hypothesis, exiting non-zero if anything required at this
  stage fails. `--at=start` or `--at=reactant` selects which criteria are
  binding, because a criterion describing the reactant state is not a defect in
  a structure that has yet to be equilibrated into one.

  `--hypothesis=` is **not optional** when more than one is defined, and
  omitting it is refused rather than defaulted. A criterion is only ever a
  criterion *within* a mechanism: "the closest substrate oxygen must be the
  β-lactam carbonyl" is a fact about the carbonyl-on-metal mechanism, which is
  the thing under test, and a gate that rejects a rival hypothesis by assuming
  the favoured one is not a gate (D127).

  **A coordination criterion returns a count, and a count is not a
  composition.** Four donors inside the cutoff can be the wrong four. Where the
  identity of the sphere is the claim, read it off `shell`'s occupancies rather
  than off `check`'s number (D129, D130).

## Why the selector check exists

`sections.py` resolves every `chain:resid:atom` in a section against that
section's own start structure and refuses to write if any is missing. Two days
of compute were lost to a 2D grid whose second coordinate began 4.35 Å apart —
which drags a group across space rather than transferring a proton. Check the
coordinate against the structure it will actually start from.

`restraints` is in that list as of D128, and it was the one place missing: a
section whose only selectors were its restraint pairs printed `0 selectors
checked` and was written anyway, while `md.md_restraints` resolves those same
pairs hours into the run.

## Email notifications

`run.slurm` carries `#SBATCH --mail-type=END,FAIL` but deliberately **no
`--mail-user`**: this repository is public, so an address committed here would
be scraped. Pass one at submit time, which overrides the directives:

```bash
sbatch --mail-user=you@example.com ... run.slurm ...
```

For a job already queued or running: `scontrol update JobId=<id>
MailUser=you@example.com MailType=END,FAIL`.
