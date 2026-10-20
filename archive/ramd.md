# τRAMD in mcube: design notes

`mcube ramd` estimates relative residence times with the τRAMD protocol. This note records how the implementation follows the published method and where it had to choose. The config keys are in the README.

## Sources

- **RAMD:** Lüdemann, Lounnas & Wade, *J. Mol. Biol.* 2000.
- **τRAMD protocol and analysis:** Kokh et al., *J. Chem. Theory Comput.* 2018. Reference script: `tauRAMD.py` (HITS-MCM/tauRAMD).
- **GROMACS implementation and recommended parameters:** Kokh et al., *J. Chem. Phys.* 2020 (arXiv:2006.11066).
- **seekrcentral/openmm_ramd:** read as a reference implementation only. It is early-stage and does not compute residence times.

## What follows the published protocol

| Item | Published | mcube |
| --- | --- | --- |
| Force | on the ligand COM, 14 kcal/mol/Å to start | `force: 14` on the binder COM (mass-weighted); fixed for a whole comparison set |
| Direction check | every 100 fs | `eval_interval: 0.1` ps |
| Redraw threshold | ligand COM displacement **relative to the protein COM** < 0.025 Å | `rmin: 0.025`, measured relative to the receptor COM |
| Stop | ligand–protein COM distance > protein extent + 5–10 Å | `max_distance`, default receptor radius + binder radius + 10 Å |
| Thermostat | Langevin, 1 ps relaxation (NAMD); the paper notes it strongly affects dissociation times | LangevinMiddle, `friction: 1.0` /ps |
| Ensemble | NPT (Langevin barostat in NAMD) | Monte Carlo barostat at 1 bar (explicit solvent) |
| Replicas | ≥ 4 starting snapshots from independent equilibrations, 15–30 trajectories each | `replicas: 4`, `trajectories: 15`; per-replica unbiased equilibration with new velocities; `replica_pdbs` for independent snapshots |
| τ per replica | 50% point of the dissociation CDF, bootstrapped | Mean of the medians of random 80% subsets (as `tauRAMD.py`), 5000 rounds |
| Overall τ | mean ± SD over replicas | Same (population SD, as `tauRAMD.py`) |
| KS check | vs Poisson CDF with τ, 0.1–0.2 typical | Same model, exact ECDF instead of `tauRAMD.py`'s binned CDF |

## Choices to revisit

1. **Binder radius in `max_distance`.** The paper's rule (protein extent + 5–10 Å) assumes a small ligand. For protein binders, the binder's own radius is added so that contact is really lost. Small ligands are barely affected.
2. **Centre-of-mass motion removal.** The system keeps OpenMM's CMMotionRemover (`remove_cm_motion: true`). The RAMD force gives the system net momentum, which the remover cancels with a uniform counter-acceleration. That is negligible for a solvated protein, but **untested on a calibration set**. Compare `remove_cm_motion: false` on a few complexes.
3. **No opposite force on the receptor.** This follows NAMD and the paper's description. GROMACS-RAMD is built on the pull code, which may act on both groups; that hasn't been confirmed.
4. **Censoring.** Trajectories that hit `max_time` are kept at `max_time`. `tauRAMD.py` instead pads replicas to 15 values with a fixed 55 ns and skips replicas with ≤ 8 events. mcube flags a replica whose τ is only a lower bound (half or more censored).
5. **Monitoring.** The COM distance and the COM vector are read as CustomCVForce collective variables, so no positions are transferred at each check. The COM vector uses raw (unwrapped) coordinates; a jump of more than half a box is treated as a periodic wrap, not as slow motion.
6. **Force calibration.** Scan the force once on a calibration set: the shortest τ should exceed ~100 ps (mcube warns below that). Then freeze it. Protein–protein correlations were reported to degrade at higher forces.

## Tests

`tests/test_ramd.py` (needs only OpenMM, numpy and pytest) checks:
- the net binder force equals F·d;
- the monitored CVs match COMs computed from the positions;
- a strong force dissociates and no force does not;
- trajectories restart from the same state;
- the CMMotionRemover option works;
- the analysis recovers ln2·τ for exponential stop times;
- a censored majority is flagged.
