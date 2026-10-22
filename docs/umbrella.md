# Umbrella sampling with mcube: a practical guide

This guide explains how to get a free-energy profile (a **PMF**, potential of mean force) for an
ion moving through a membrane channel. It assumes you can already run `mcube preprocess` and
`mcube md`. The complete list of config keys is in the README (`### umbrella.inputs`); a working
config is [`config_examples/umbrella_config.yaml`](../config_examples/umbrella_config.yaml).

## What it computes, in one paragraph

An ion rarely visits the high-energy parts of a channel on its own, so plain MD cannot tell you
how hard it is to get through. Umbrella sampling splits the path into **windows**. In each window
a spring holds the ion near one height `z0` above or below the channel's centre, and the ion is
simulated there. Each window samples a little stretch of the path well. **WHAM** then removes the
effect of the springs and joins the windows into one profile, F(z), in kcal/mol. The highest
point of F(z) relative to bulk water is the barrier to permeation.

`z` is always the ion's height relative to the centre of mass of the atoms you name as the
`reference` (normally the channel's Cα atoms), along the membrane normal (the z axis).

## The two levels of theory

| leg | model | where | cost | config section |
|---|---|---|---|---|
| MM | classical force field | the whole path, bulk to bulk | GPU, ns per window | `umbrella` |
| QM/MM | GFN2-xTB for the ion and its neighbours, MM for the rest | only the narrowest part | CPU, ~20 ps per window | `umbrella_qmmm` |

MM can sample long enough to capture the ion losing and regaining water molecules, which takes
nanoseconds. QM/MM describes the chemistry of the ion binding to the protein properly, but can
only afford picoseconds. So MM gives the whole profile and QM/MM replaces it where the chemistry
matters. `stitch` joins the two, and `wham --reweight dft` can then correct the xTB part to DFT.

## Before you start

1. **An equilibrated system with the ion in it.** The channel sits in a bilayer with its axis
   along z, and the ion has replaced one water in bulk solution. Run `preprocess` (with the ion
   listed under `metals`) and a short `md` on it. `md`'s last frame is where the pull starts.
2. **Pick the metal parameters once, in `md.inputs.metal_parameters`.** Every later stage uses the
   same set (see *Metal parameters* below).
3. **Choose the reference atoms.** They must stay together, for example `"chain A and name CA"`
   for a barrel. mdtraj selection syntax; PDB residue numbers are `resSeq`, not `resid`.
4. **Choose the range.** `z_start`/`z_end` must reach bulk water on both sides (for a ~25 Å thick
   membrane, ±35 Å is typical) with `spacing: 0.5` Å.

## Step by step

Each command is its own cluster job. Run them from one working directory; each window makes its
own subdirectory (`umbrella_037/`, `umbrella_qmmm_007/`, ...) and records itself in the config.

**1. Seed the windows (one GPU job).**

```bash
mcube umbrella config.yaml pull
```

A spring is dragged across the whole range at `pull_speed` (Å/ns), and a structure is saved as
it passes each window centre. Check `pull/pull_cv.csv`: the column `z_A` should follow `z0_A`
with a small lag. If it lags by several Å, the ion is stuck; use a slower pull or a stiffer `pull_k`.

**2. Run every MM window (one GPU job per window).**

```bash
mcube umbrella config.yaml window --window 0
mcube umbrella config.yaml window --window 1
...
```

Windows are numbered from 0 at `z_start`. Each one minimises, equilibrates for `equilibrate` ps
with the ion's first hydration shell held, then lets go and samples for `time` ps.

**3. Build the MM profile.**

```bash
mcube umbrella config.yaml wham
```

This writes `umbrella_pmf.csv` (columns: z, F, error) and records it under `umbrella.results.pmf`.
Look at the warning line: if neighbouring windows hardly overlap, add a window between them or
lower `k` there, and run `wham` again.

**4. Run the QM/MM windows over the constriction (one CPU job per window).**

First decide `qm_residues` in `umbrella_qmmm`: the protein residues that touch the ion anywhere
in the QM/MM range, read from the MM windows there. Then:

```bash
mcube umbrella config.yaml window --window 0 --section umbrella_qmmm
...
mcube umbrella config.yaml wham --section umbrella_qmmm
```

Each QM/MM window starts from the finished MM window at the same `z0`.

**5. Join the two profiles.**

At the first and last QM/MM window centres, score frames from both the MM window and the QM/MM
window at both levels, then stitch. The MM window at height z is number
`(z - z_start) / spacing` of the `umbrella` section (z = -8 Å with `z_start: -35` and
`spacing: 0.5` is window 54):

```bash
mcube umbrella config.yaml energies --window <first MM index at that z> --level mm
mcube umbrella config.yaml energies --window <same>                     --level xtb
mcube umbrella config.yaml energies --window 0 --section umbrella_qmmm  --level mm
mcube umbrella config.yaml energies --window 0 --section umbrella_qmmm  --level xtb
# ... the same four for the last QM/MM centre ...
mcube umbrella config.yaml stitch
```

This writes `stitched_umbrella_qmmm.csv`, the combined profile.

**6. (Optional) Correct xTB to DFT.**

```bash
mcube umbrella config.yaml energies --window N --section umbrella_qmmm --level dft --stride 5
mcube umbrella config.yaml energies --window N --section umbrella_qmmm --level xtb --stride 5
mcube umbrella config.yaml wham --section umbrella_qmmm --reweight dft
```

Use the same `--stride` for both levels: frames are matched one to one. DFT is expensive, so
start with a large stride.

## How to tell whether a result is trustworthy

| check | where | what you want |
|---|---|---|
| window overlap | `wham` warning, `pmf.min_overlap` | above 0.05 between every pair of neighbours |
| statistical error | `err_kcal` column of the PMF | well below the barrier you want to quote |
| convergence | extend the windows (larger `time`, rerun) and run `wham` again | the profile stops moving, within the error |
| bulk is flat | both ends of the PMF | F(z) levels off; if not, extend the range |
| stitch closure | `umbrella_qmmm.results.stitch.closure_kcal` | within ~1 kcal/mol; a larger value shows as a step in the stitched profile |
| QM region intact | `shell_changes` in each QM/MM window entry | empty, or only second-shell waters |

A large stitch closure means one of two things: the QM/MM windows are too short, or the
correction from MM to QM/MM changes beyond the ends of the QM/MM range. Lengthen the QM/MM
windows first; if the closure persists, widen the QM/MM range.

## Metal parameters

A classical force field has no chemical bonds to a metal ion, so the ion's size parameters decide
how closely oxygens pack around it. `metal_parameters` chooses those parameters:

| value | what it is | use it for |
|---|---|---|
| `iod` (default) | Li/Merz 12-6 parameters fitted to ion–oxygen distances | any metal; the safe choice |
| `1264` | Li/Song/Merz 12-6-4: adds an attraction to every atom that depends on that atom's polarisability; better for highly charged ions | +3/+4 ions such as La³⁺ and Lu³⁺ (not yet +2 ions such as Ca²⁺) |
| `none` | whatever the force-field file contains | a control for Zn²⁺ only. **Never for lanthanides:** OpenMM's stock La³⁺/Lu³⁺ parameters are half of the 12-6-4 model |

Set it once in `md.inputs`. Since this change, `qmmm`, `scan`, `neb`, `ts`, `sp` and `umbrella`
all use the set `md` used, unless their own section says otherwise. Before this, the QM/MM stages
silently used the force field's own metal parameters, so QM/MM results from before and after
the change are not directly comparable.

Under `1264`, every atom type needs a polarisability. Protein, water, ion and Lipid17 atoms are
covered. A ligand's atom types are not: give them as
`polarizabilities: {<atom class>: <Å³>}` in `md.inputs`, or mcube stops and names the missing types.

## Common mistakes

- **Using `resid` in the reference selection.** mdtraj's `resid` counts from 0; use `resSeq` for
  the residue numbers you see in the PDB file.
- **A reference group that spans the periodic box.** mcube refuses to start. Give a structure in
  which the channel is whole.
- **Changing the QM region partway through the QM/MM leg.** `qm_waters` and `qm_residues` must be
  the same in every QM/MM window, or the stitch compares different models.
- **Different `--stride` for `dft` and `xtb` energies.** `wham --reweight` refuses mismatched frames.
- **Reading a profile before `wham` stops warning about overlap.** Gaps make F(z) unreliable there.

## Where things are written

| file | written by | contents |
|---|---|---|
| `pull/seed_NNN.pdb`, `pull/pull_cv.csv` | `pull` | the starting structure of each window; the pull trace |
| `<section>_NNN/cv.csv` | `window` | the ion's z (and sideways distance r) over time |
| `<section>_NNN/window.dcd`, `lastframe.pdb` | `window` | saved frames; the final structure |
| `<section>_NNN/energies_<level>.csv` | `energies` | the energy of each saved frame at that level |
| `<section>_pmf.csv` | `wham` | z, F (kcal/mol), bootstrap error |
| `stitched_<section>.csv` | `stitch` | the combined MM + QM/MM profile |

Everything is also recorded in the config under `<section>.results`, so a later step finds what
an earlier one produced without being told the paths.

## Further reading

- Kästner, *WIREs Comput. Mol. Sci.* 2011, 1, 932 — an introduction to umbrella sampling.
- Allen, Andersen & Roux, *Biophys. J.* 2006, 90, 3447 — PMFs for ions in channels: setup and convergence.
- Heimdal & Ryde, *PCCP* 2012, 14, 12592 — why the QM/MM leg is sampled directly rather than
  computed by perturbing the MM leg.
- Li et al., *J. Chem. Theory Comput.* 2018, 14, 5583 — reweighting a cheap QM/MM profile to a
  higher level.
- Li & Merz, *Chem. Rev.* 2017, 117, 1564 — metal ion models, including 12-6-4.
