# Choosing the QM level: xTB or DFT

Every gradient-driven stage (`qmmm` optimisation, `scan`, `neb`, `ts`, and the
frequencies inside `ts`) builds a QM/MM theory. By default the QM part is
GFN2-xTB. Setting `qm_engine: dft` runs those same stages with PySCF instead.

The `singlepoint` at the end of the `qmmm` stage is always DFT and is unaffected
by `qm_engine`; it has its own `sp_*` keys.

## Which to use

| | GFN2-xTB (`qm_engine: xtb`) | DFT (`qm_engine: dft`) |
| --- | --- | --- |
| Cost per gradient | seconds | minutes to hours |
| Good for | exploring: MD, scans, NEB, locating approximate stationary points | refining a structure an xTB path already found, and the numbers you publish |
| Not good for | quantitative barriers | anything that needs many gradients |

The usual pattern is both, in sequence:

1. `qmmm` with the default engine for QM/MM MD and a first minimisation.
2. `scan` / `neb` at xTB to find the path and the approximate transition state.
3. `qmmm` or `ts` with `qm_engine: dft` to re-optimise each stationary point,
   get frequencies and run the IRC.
4. Single points at a larger basis on those geometries for the final energies.

Do not put `qm_engine: dft` on a scan or NEB unless you have costed it: a
25-point scan is at least 25 optimisations, each of many gradients.

## How the level is specified

The DFT engine reads `dft_*` keys, and **each one falls back to its `sp_*`
counterpart**. So if a config already defines the single-point level, the
gradients use the same level with no extra keys:

```yaml
qmmm:
  inputs:
    sp_functional: wb97x
    sp_basis: def2-svp
    sp_gridlevel: 3
    sp_densityfit: true

# Gradients at wB97X/def2-SVP too, because dft_* falls back to sp_*
qmmm_refine:
  inputs:
    qm_engine: dft
```

To use a cheaper level for the geometry than for the energy — the common choice,
because geometries converge at a lower level than energies — set both:

```yaml
qmmm_refine:
  inputs:
    qm_engine: dft
    dft_functional: r2scan          # meta-GGA: no exact exchange, much cheaper
    dft_scf_dispersion: D4          # dispersion inside the SCF, so it is in the gradient
    dft_basis: def2-svp
    sp_functional: wb97x            # the reported energy, on the r2SCAN geometry
    sp_basis: def2-tzvpd
    sp_dispersion: d3bj
```

Note the two different dispersion keys, which do different things:

- `dft_scf_dispersion` is ASH's own `dispersion` keyword (`D3`, `D4`, `VV10`, …).
  It is part of the SCF, so it enters the **gradient** and therefore affects the
  optimised geometry. Use this when optimising.
- `sp_dispersion` (`d3bj` / `d3zero`) is mcube's own correction, computed with
  `simple-dftd3` **after** a single point and reported separately as
  `sp_disp_energy` / `sp_energy_disp`. Use this when you want the dispersion
  contribution as its own number.

Setting both is fine — they apply to different steps — but do not set
`dft_scf_dispersion` *and* choose a functional that already includes dispersion,
or you will count it twice.

## Where to put the keys

`qm_engine` and the `dft_*` / `sp_*` keys can go on `qmmm.inputs`, where every
stage picks them up, or on an individual `scan` / `neb` / `ts` section, which
overrides them for that stage only. That matters because those three stages take
the QM theory from `qmmm.inputs`: before this was wired up, a `qm_engine: dft` on
a `ts_*` section was silently ignored and the stage quietly ran at xTB.

Everything else about the QM/MM setup — charge, multiplicity, force fields, the
QM region — always comes from `qmmm.inputs`, so a stage section cannot change the
chemistry by accident. (The QM region itself is the exception: `scan`, `neb` and
`ts` accept `qm_residues` / `extra_qm_residues` per section.)

So a typical mechanism config keeps the exploring level at the top and names the
refinement level only on the sections that refine:

```yaml
qmmm:
  inputs:
    sp_functional: wb97x        # single points, and the fallback for dft_*
    sp_basis: def2-svp

scan_attack:                    # xTB, because qm_engine is not set anywhere
  inputs: {atoms: [...], end: 1.45, points: 25}

ts_step1:                       # this stage, and only this stage, at DFT/MM
  inputs:
    qm_engine: dft
    dft_functional: r2scan
    sp_basis: def2-tzvpd
```

## Frequencies: use a partial Hessian

The `ts` stage computes a numerical Hessian, which costs `3N × 2` gradients for
`N` atoms. At xTB that is routine. At DFT it is the most expensive thing in the
whole protocol:

| QM atoms in the Hessian | Gradients | At 1 min each | At 30 min each |
| --- | --- | --- | --- |
| 120 | 720 | 12 h | 15 days |
| 25 | 150 | 2.5 h | 3 days |
| 15 | 90 | 1.5 h | ~2 days |

So with `qm_engine: dft`, always restrict the Hessian to the atoms that actually
move in the reaction:

```yaml
ts:
  inputs:
    qm_engine: dft
    partial_hessian_atoms: ["C:2:O01", "C:1:C4x", "C:1:O2x", "C:1:N1x", "A:226:HZ1"]
```

The imaginary mode is still the reaction mode; what you lose is the other
frequencies, which you were not going to quote anyway. If you need a full
Hessian, do it at xTB and report it as such.

## GPU

PySCF can run on a GPU through `gpu4pyscf`, which is a large speed-up for a QM
region of ~100 atoms:

```yaml
qmmm_refine:
  inputs:
    qm_engine: dft
    dft_platform: GPU
    dft_gpu_pcgrad: true     # also put the MM point-charge gradient on the device
```

`dft_gpu_pcgrad` matters for QM/MM specifically: with electrostatic embedding the
MM environment enters as point charges, and for a solvated system that gradient
is a large share of the total. It is off by default because it is off in ASH.

Installing `gpu4pyscf` needs care, because the current `cupy` requires numpy ≥ 2
and will happily upgrade numpy out from under OpenMM and mdtraj, which are
compiled against numpy 1.x. Check first and pin what must not move:

```bash
pip install --dry-run gpu4pyscf-cuda12x
```

Read the `Would install` line. If it lists `numpy`, constrain it:

```bash
pip install gpu4pyscf-cuda12x "cupy-cuda12x==13.3.0" "numpy==1.26.4"
```

Then confirm the environment still works:

```bash
python -c "import numpy, openmm, mdtraj, pyscf; print(numpy.__version__, openmm.__version__)"
```

Use the `cuda11x` package instead if that is what the machine's driver supports.

## Cores and memory

`dft_cores` sets the thread count for the DFT gradients, defaulting to the
stage's `cores`. `dft_memory` is PySCF's `max_memory` in MB and defaults, like
`sp_memory`, to 80% of `SLURM_MEM_PER_NODE` so density-fitting integrals stay in
RAM instead of going to disk; ASH's own default of 3100 MB is usually too small
and makes a large QM region much slower.

Thread binding (see `running-on-bunya.md`) only ever applies to xTB, which runs
as a separate process. The DFT engine runs in-process, so it is deliberately left
unbound — binding the whole process would pin PySCF's BLAS threads to one core.

## Waters in the QM region

While optimising, MM waters in the active region are held rigid, but **waters in
the QM region are left free**, so a QM water can donate or accept a proton. If
you are studying a nucleophilic water or a proton relay, put those waters in the
QM region (`extra_qm_residues`) and they will be able to move their hydrogens.
`rigid_water: false` in the `ts` stage frees the MM ones too.
