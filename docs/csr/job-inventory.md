# Every CSR-1 job that has run, and what came out of it

One row per Slurm allocation, from `sacct` rather than from memory, with the
decision-log entries that read it. **Regenerate it the same way:**

```
sacct -S 2026-08-01 -X -o jobid,jobname%16,state,elapsed,workdir%44 \
  | grep CSR
```

`-X` gives one line per allocation rather than per step. The `D` column is
derived by matching job numbers against `docs/csr-decision-log.md`, so an
empty cell means the log never cites that job - which for a setup failure is
usually right and for a completed run is a gap worth closing.

## Where the run directories are, since 2026-11-03

All thirteen now live under `/scratch/user/uquadhik/CSR/`, and **nothing** is
left at the top of scratch: `CSR/CSR_wt`, `CSR/CSR_crystal_v2`, `CSR/CSR_nac`,
`CSR/CSR_mich` and so on.

Twelve moved at 20:50. **`CSR_mich` followed at 22:25**, once its MD (31721717)
and then its postprocess (31736226) had both finished - Slurm records a
`WorkDir` at submit time, so moving a running job's directory breaks it, and
that gate was checked before each attempt rather than assumed from an expected
finish time.

Paths written in the entries below, and in the decision log, are the paths as
they were when the work ran. They have not been rewritten - a record of what was
done should say what was done - so prefix a `CSR_*/...` path with `CSR/` when
following it, `CSR_mich` aside.

What *was* rewritten is the 508 absolute `/scratch/user/uquadhik/CSR_*` paths
inside the 13 moved `config.yaml` and `.slurm` files, which are read by the
pipeline rather than by a person, so a stale one is a silent failure. Each of
those files has a `.premove` copy beside it. `.out` logs were deliberately left
alone, for the same reason the entries below were. `site.yaml`'s `aspout`
already used `../CSR_nac/...` and so needed nothing: the runs moved together, so
a sibling is still a sibling.

A COMPLETED state means Slurm's accounting was happy, **not** that the job
produced anything: job 31186038 is in this table as COMPLETED and computed
nothing at all (¶153). That is why `mcube.cli.confirm_produced` now exits
non-zero when a stage's result carries no product key.

## Runs with something to show

| job | name | directory | state | elapsed | D | what it was / what came out |
| --- | --- | --- | --- | --- | --- | --- |
| 29835095 | `mcube_QMMM` | CSR_hydrolysed | COMPLETED | 00:01:41 |  |  |
| 29836509 | `mcube_QMMM` | CSR_ligand | CANCELLED | 00:00:00 |  |  |
| 29836510 | `mcube_QMMM` | CSR_hydrolysed | COMPLETED | 00:04:46 |  |  |
| 29837442 | `csr_qmmm_opt` | CSR_ligand | CANCELLED | 00:00:00 |  |  |
| 29837804 | `csr_qmmm_opt` | CSR_ligand | CANCELLED | 00:34:31 |  |  |
| 29840911 | `csr_md_zw` | CSR_ligand | CANCELLED | 00:31:34 |  |  |
| 29848762 | `csr_md_restr` | CSR_crystal | COMPLETED | 04:33:35 |  | first restrained MD of the crystal system |
| 29893675_0 | `csr_qmmm_probe` | CSR_crystal | COMPLETED | 00:24:26 |  | three QM/MM probes of the crystal site |
| 29893675_1 | `csr_qmmm_probe` | CSR_crystal | COMPLETED | 00:27:46 |  | three QM/MM probes of the crystal site |
| 29893675_2 | `csr_qmmm_probe` | CSR_crystal | COMPLETED | 02:09:24 |  | three QM/MM probes of the crystal site |
| 29919076 | `csr_scan_f400` | CSR_crystal | CANCELLED | 00:11:34 |  |  |
| 29920397 | `csr_v2_mm_qmmm` | CSR_crystal_v2 | CANCELLED | 04:45:01 |  |  |
| 29937049 | `csr_v2_qmmm` | CSR_crystal_v2 | CANCELLED | 00:13:16 |  |  |
| 29937358 | `xtb_bench` | CSR_crystal_v2/xtb_bench | COMPLETED | 00:01:47 |  |  |
| 29937641 | `csr_v2_qmmm` | CSR_crystal_v2 | FAILED | 18:30:02 |  | died at the first optimisation step on ASH's missing `MMenergy` after OpenMM_MD - the reason `qmmm` now rebuilds the theory |
| 29981678_0 | `csr_protdft` | CSR_crystal_v2 | CANCELLED | 03:58:43 |  |  |
| 29981678_1 | `csr_protdft` | CSR_crystal_v2 | CANCELLED | 03:58:43 |  |  |
| 29981678_2 | `csr_protdft` | CSR_crystal_v2 | CANCELLED | 03:58:42 |  |  |
| 29992285_0 | `csr_protdft2` | CSR_crystal_v2 | CANCELLED | 00:11:48 |  |  |
| 29992285_1 | `csr_protdft2` | CSR_crystal_v2 | CANCELLED | 00:11:48 |  |  |
| 29992285_2 | `csr_protdft2` | CSR_crystal_v2 | CANCELLED | 00:11:47 |  |  |
| 29996196_0 | `csr_protdft2` | CSR_crystal_v2 | CANCELLED | 00:24:41 |  |  |
| 29996196_1 | `csr_protdft2` | CSR_crystal_v2 | CANCELLED | 00:24:41 |  |  |
| 29996196_2 | `csr_protdft2` | CSR_crystal_v2 | CANCELLED | 00:24:40 |  |  |
| 29997262_0 | `csr_protdft2` | CSR_crystal_v2 | CANCELLED | 00:07:06 |  |  |
| 29997262_1 | `csr_protdft2` | CSR_crystal_v2 | CANCELLED | 00:07:06 |  |  |
| 29997262_2 | `csr_protdft2` | CSR_crystal_v2 | CANCELLED | 00:07:06 |  |  |
| 29997628_0 | `csr_protdft2` | CSR_crystal_v2 | FAILED | 10:58:18 |  | DFT single points; task 0 failed, 1 and 2 completed |
| 29997628_1 | `csr_protdft2` | CSR_crystal_v2 | COMPLETED | 09:52:04 |  | DFT single points; task 0 failed, 1 and 2 completed |
| 29997628_2 | `csr_protdft2` | CSR_crystal_v2 | COMPLETED | 10:36:31 |  | DFT single points; task 0 failed, 1 and 2 completed |
| 30022363_0 | `csr_scans` | CSR_crystal_v2 | OUT_OF_MEMORY | 04:44:25 |  | 12-task scan array, ten of them killed for memory |
| 30022363_1 | `csr_scans` | CSR_crystal_v2 | COMPLETED | 01:24:21 |  | 12-task scan array, ten of them killed for memory |
| 30022363_2 | `csr_scans` | CSR_crystal_v2 | OUT_OF_MEMORY | 03:15:31 |  | 12-task scan array, ten of them killed for memory |
| 30022363_3 | `csr_scans` | CSR_crystal_v2 | CANCELLED | 11:08:24 |  | 12-task scan array, ten of them killed for memory |
| 30022363_4 | `csr_scans` | CSR_crystal_v2 | CANCELLED | 09:37:15 |  | 12-task scan array, ten of them killed for memory |
| 30022363_5 | `csr_scans` | CSR_crystal_v2 | OUT_OF_MEMORY | 08:20:55 |  | 12-task scan array, ten of them killed for memory |
| 30022363_6 | `csr_scans` | CSR_crystal_v2 | OUT_OF_MEMORY | 06:24:54 |  | 12-task scan array, ten of them killed for memory |
| 30022363_7 | `csr_scans` | CSR_crystal_v2 | OUT_OF_MEMORY | 07:32:22 |  | 12-task scan array, ten of them killed for memory |
| 30022363_8 | `csr_scans` | CSR_crystal_v2 | OUT_OF_MEMORY | 02:59:04 |  | 12-task scan array, ten of them killed for memory |
| 30022363_9 | `csr_scans` | CSR_crystal_v2 | OUT_OF_MEMORY | 04:25:46 |  | 12-task scan array, ten of them killed for memory |
| 30022363_10 | `csr_scans` | CSR_crystal_v2 | OUT_OF_MEMORY | 01:51:19 |  | 12-task scan array, ten of them killed for memory |
| 30022363_11 | `csr_scans` | CSR_crystal_v2 | OUT_OF_MEMORY | 03:40:38 |  | 12-task scan array, ten of them killed for memory |
| 30047537_0 | `csr_react` | CSR_crystal_v2 | CANCELLED | 00:07:21 |  |  |
| 30047537_1 | `csr_react` | CSR_crystal_v2 | CANCELLED | 00:07:21 |  |  |
| 30047537_2 | `csr_react` | CSR_crystal_v2 | CANCELLED | 00:07:21 |  |  |
| 30047537_3 | `csr_react` | CSR_crystal_v2 | CANCELLED | 00:07:21 |  |  |
| 30047537_4 | `csr_react` | CSR_crystal_v2 | CANCELLED | 00:07:21 |  |  |
| 30047537_5 | `csr_react` | CSR_crystal_v2 | CANCELLED | 00:07:21 |  |  |
| 30047537_6 | `csr_react` | CSR_crystal_v2 | CANCELLED | 00:07:21 |  |  |
| 30047537_7 | `csr_react` | CSR_crystal_v2 | CANCELLED | 00:07:21 |  |  |
| 30047537_8 | `csr_react` | CSR_crystal_v2 | CANCELLED | 00:07:21 |  |  |
| 30048252_0 | `csr_react` | CSR_crystal_v2 | OUT_OF_MEMORY | 01:15:41 |  | reactant-state array, four completed |
| 30048252_1 | `csr_react` | CSR_crystal_v2 | COMPLETED | 00:12:19 |  | reactant-state array, four completed |
| 30048252_2 | `csr_react` | CSR_crystal_v2 | COMPLETED | 00:09:12 |  | reactant-state array, four completed |
| 30048252_3 | `csr_react` | CSR_crystal_v2 | COMPLETED | 00:08:29 |  | reactant-state array, four completed |
| 30048252_4 | `csr_react` | CSR_crystal_v2 | CANCELLED | 01:22:58 |  | reactant-state array, four completed |
| 30048252_5 | `csr_react` | CSR_crystal_v2 | COMPLETED | 00:28:12 |  | reactant-state array, four completed |
| 30048252_6 | `csr_react` | CSR_crystal_v2 | CANCELLED | 01:22:58 |  | reactant-state array, four completed |
| 30048252_7 | `csr_react` | CSR_crystal_v2 | CANCELLED | 01:22:58 |  | reactant-state array, four completed |
| 30048252_8 | `csr_react` | CSR_crystal_v2 | COMPLETED | 00:56:06 |  | reactant-state array, four completed |
| 30048856 | `scan2d_H1_f00497` | CSR_crystal_v2 | CANCELLED | 00:13:51 |  |  |
| 30053828 | `csr_cleanR` | CSR_crystal_v2 | COMPLETED | 09:51:17 | ¶47, ¶48, ¶60, ¶70, ¶112 | `cleanR` - the clean reactant QM/MM optimisation |
| 30054626 | `csr_gpubench` | CSR_crystal_v2 | COMPLETED | 00:01:14 | ¶47 | GPU benchmark set (four runs) |
| 30055503 | `csr_gpubench` | CSR_crystal_v2 | COMPLETED | 00:01:16 | ¶47 |  |
| 30056056 | `csr_gpubench` | CSR_crystal_v2 | COMPLETED | 00:01:16 | ¶47 |  |
| 30056429 | `csr_gpubench` | CSR_crystal_v2 | COMPLETED | 00:01:15 | ¶47 |  |
| 30085418 | `csr_md2` | CSR_crystal_v2 | COMPLETED | 09:11:36 | ¶50, ¶51, ¶60 | `md2` - production MD behind ¶50/¶51 |
| 30119957 | `csr_zncheck` | CSR_crystal_v2 | COMPLETED | 1-01:10:50 | ¶53, ¶54, ¶56, ¶57, ¶58, ¶61 | Zn-Asp92 walk plus DFT single points, 25 h; ¶53-¶61 |
| 30122134 | `csr_dftprobe` | CSR_crystal_v2 | FAILED | 00:01:08 |  |  |
| 30123101 | `csr_dftprobe` | CSR_crystal_v2 | CANCELLED | 01:59:11 | ¶57 |  |
| 30200387 | `csr_coscan` | CSR_crystal_v2 | COMPLETED | 00:49:49 | ¶58, ¶59, ¶61 | carbonyl scan (`coscan`) |
| 30205872 | `v3_md` | CSR_crystal_v3 | CANCELLED | 00:06:27 | ¶62, ¶63 |  |
| 30207244 | `spco0` | CSR_crystal_v2 | COMPLETED | 00:26:17 | ¶63, ¶70 | wB97X-D3/def2-SVP single points on coscan points 0, 4, 8 |
| 30207245 | `spco4` | CSR_crystal_v2 | COMPLETED | 00:25:33 | ¶63, ¶70 |  |
| 30207246 | `spco8` | CSR_crystal_v2 | COMPLETED | 00:25:54 | ¶63, ¶70 |  |
| 30226271 | `umd` | CSR_unrestrained | COMPLETED | 20:44:49 | ¶69, ¶71, ¶72, ¶73 | 100 ns unrestrained MD - the run that showed the shell drifts |
| 30226272 | `uqmmm` | CSR_unrestrained | CANCELLED | 01:52:05 | ¶69, ¶71, ¶72 |  |
| 30302149 | `g93md` | CSR_unrestrained2 | CANCELLED | 01:47:37 |  |  |
| 30312014 | `r6prep` | CSR_unrestrained2 | COMPLETED | 00:01:23 |  |  |
| 30312015 | `r6md` | CSR_unrestrained2 | FAILED | 00:00:10 |  |  |
| 30313764 | `r6prep` | CSR_unrestrained2 | COMPLETED | 00:02:13 |  |  |
| 30313765 | `r6md` | CSR_unrestrained2 | CANCELLED | 00:01:52 |  |  |
| 30314534 | `g6prep` | CSR_unrestrained2 | COMPLETED | 00:01:28 |  |  |
| 30314535 | `g6md` | CSR_unrestrained2 | FAILED | 00:00:13 |  |  |
| 30348138 | `f6prep` | CSR_unrestrained2 | COMPLETED | 00:01:28 |  |  |
| 30348139 | `f6md` | CSR_unrestrained2 | FAILED | 00:00:18 |  |  |
| 30349570 | `h6prep` | CSR_unrestrained2 | COMPLETED | 00:01:46 |  |  |
| 30349571 | `h6md` | CSR_unrestrained2 | FAILED | 01:07:40 |  |  |
| 30354772 | `p6md` | CSR_unrestrained2 | COMPLETED | 02:35:00 |  | pose-6 MD and QM/MM |
| 30354773 | `p6qm` | CSR_unrestrained2 | COMPLETED | 01:44:14 |  |  |
| 30394832 | `gln93` | CSR_gln93 | COMPLETED | 00:08:58 | ¶93 | Gln93 system build |
| 30409808 | `zngln` | CSR_gln93 | OUT_OF_MEMORY | 00:59:21 | ¶93, ¶97, ¶105 | `zngln` scan, killed for memory - why scan memory is now 64G |
| 30416229 | `nacprep` | CSR_nac | COMPLETED | 00:01:44 | ¶98 | `nac` build |
| 30416230 | `nacmd` | CSR_nac | CANCELLED | 01:10:38 | ¶98, ¶99 |  |
| 30418017 | `nacgpu` | CSR_nac | COMPLETED | 00:04:20 | ¶100 | GPU MD rate check (0.52 ns/min) |
| 30421166 | `naclong` | CSR_nac | COMPLETED | 06:19:45 | ¶101, ¶102, ¶105, ¶128, ¶134 | `naclong` - 500 ps NAC MD; its last frame fed everything after |
| 30466101 | `nacqm` | CSR_nac | CANCELLED | 00:02:25 | ¶103 |  |
| 30467248 | `nacqm` | CSR_nac | CANCELLED | 1-09:44:28 | ¶103, ¶104, ¶105, ¶151 | `nacqm` - 33 h of QM/MM MD, cancelled; source of ¶151's soft 1.7 ps/h figure |
| 30580181 | `nacopt` | CSR_nac | CANCELLED | 00:00:24 |  |  |
| 30580259 | `nacopt` | CSR_nac | FAILED | 00:59:52 | ¶106 | `nacopt`, failed |
| 30598032 | `nacarg` | CSR_nac | COMPLETED | 23:56:12 | ¶109, ¶110, ¶112 | `nacarg` - 24 h QM/MM MD with Arg229 added to the QM region |
| 30698605 | `ohHmin` | CSR_nac | CANCELLED | 00:00:00 | ¶114, ¶115 |  |
| 30699162 | `ohHmin` | CSR_nac | COMPLETED | 00:03:16 | ¶115 | hydroxide-protonated minimisation (`ohHmin`) |
| 30703062 | `ohHmin2` | CSR_nac | FAILED | 00:00:09 | ¶116 |  |
| 30705703 | `ohHmin3` | CSR_nac | COMPLETED | 00:03:53 | ¶116, ¶118 | `ohHmin3` - the build that `ohH` ran from |
| 30721365 | `ohH` | CSR_nac | COMPLETED | 19:46:57 | ¶118, ¶119 | `ohH` - 20 h QM/MM MD of the Lys226-protonated species |
| 30923516 | `subst` | CSR_nac | CANCELLED | 02:11:54 | ¶128, ¶129, ¶130, ¶131, ¶135, ¶145 | `subst` - cancelled; its start structure carried the shell ¶120 condemned |
| 30966977 | `wt100` | CSR_wt | COMPLETED | 03:15:14 | ¶132, ¶134 | `wt100` - 100 ns, substrate left the site |
| 31016877 | `wt100b` | CSR_wt | COMPLETED | 03:12:51 | ¶136, ¶137, ¶138 | `wt100b` - 100 ns with the substrate tethered; six-coordinate zinc in 4,998 of 5,000 frames |
| 31016878 | `aspout` | CSR_wt | CANCELLED | 04:38:57 | ¶136, ¶137, ¶139, ¶140, ¶141, ¶142 | `aspout` cancelled at point 0 of 8: rigid-water constraints written as `$set` can never satisfy `convergence_cmax` |
| 31051634 | `aspout` | CSR_wt | COMPLETED | 01:16:44 | ¶139, ¶140, ¶141 | `aspout` again with `rigid_water: false` - 9/9 converged in 1 h 17 m against 0/9 in 4 h 39 m |
| 31143211 | `aspin` | CSR_wt | COMPLETED | 00:53:22 | ¶140, ¶142, ¶143 | `aspin` - the inward arm; `scan_00` (3.60 A) and `scan_08` (1.90 A) are the structures everything since has used |
| 31146622 | `rwoff` | CSR_wt | COMPLETED | 00:03:58 | ¶142 | `rwoff`/`rwon` - the one-variable test of holding the waters (+0.0027 kcal/mol) |
| 31146643 | `rwon` | CSR_wt | COMPLETED | 00:04:04 | ¶142 |  |
| 31154229 | `mich2d` | CSR_wt | CANCELLED | 00:11:17 | ¶144 | first `mich2d`, cancelled for starting from the over-coordinated structure |
| 31155441 | `mich2d` | CSR_wt | FAILED | 03:26:23 | ¶145, ¶146, ¶147, ¶148, ¶149 | the 45-point 2D grid: all 45 points converged, then the job FAILED writing its own summary (`target1_A`) |
| 31160947 | `substopt` | CSR_nac | FAILED | 00:26:44 | ¶146, ¶147, ¶152 | `substopt` - died in `xtbpcgradientgrab` on a point charge 0.58 A from a QM atom |
| 31167021 | `subst` | CSR_nac | CANCELLED | 00:00:00 | ¶147, ¶148, ¶149 |  |
| 31176144 | `spref` | CSR_wt | FAILED | 00:00:08 | ¶152 | the charge-deletion pair, dead in 8 s on a doubled section prefix (`sp_sp_qmmm_ref`) |
| 31176145 | `spr229` | CSR_wt | FAILED | 00:00:08 | ¶152 |  |
| 31185324 | `spref` | CSR_wt | COMPLETED | 03:01:51 | ¶152, ¶153 | the pair re-run: -5076.869006116424 and -5076.939666070885 Eh, Arg229 at +44.34 kcal/mol |
| 31185325 | `spr229` | CSR_wt | COMPLETED | 03:06:14 | ¶152, ¶153 |  |
| 31186037 | `asph` | CSR_wt | COMPLETED | 00:12:03 | ¶152, ¶153 | `asph` - protonated Asp92, free optimisation: Zn-OD1 1.90 -> 2.99 A, five ligands -> four |
| 31186038 | `qmrate` | CSR_wt | COMPLETED | 00:00:14 | ¶152, ¶153 | asked for 3 ps of MD, computed nothing, recorded COMPLETED |
| 31277149 | `q93` | CSR_wt | RUNNING | 00:11:09 | ¶153 | `q93` - Gln93 charge deletion |
| 31277150 | `qmrate` | CSR_wt | RUNNING | 00:11:09 | ¶153 | `qmrate` - 3 ps timing run, resubmitted |
| 31721628 | `mich_pre` | CSR_mich | COMPLETED | 00:00:51 | ¶196 | the rebuild's preprocess: `removed D:1:HOH (3 atoms)`, sphere = carbonyl 1.94, His164 2.09, OH 2.14, His88 2.22, Asp92 OD2 2.79 |
| 31721717 | `mich_eq5` | CSR_mich | submitted | - | ¶196 | 10 ns NVT at 2 fs, the three Michaelis restraints at 500 kcal/mol/A^2 |

**This table is stale between 31277150 and 31721628** - the 100 ns gate runs
(31483853 reactant, 31423830 product) and everything else from ¶186-¶195 are
missing. Regenerate it with the `sacct` command at the top before citing it as
complete.

## The `mcube_QMMM` setup-failure cluster

**40 jobs**, all named `mcube_QMMM`, all FAILED inside two minutes, in
`CSR_ligand` (33) and `CSR_hydrolysed` (7), between job 27788559 and 29835094.
None produced output and none is cited in the decision log: this is the
period when the pipeline itself was being made to run, and the failures were
read from the terminal and fixed within minutes. They are counted here
rather than listed so the table above stays readable.

Job ids: 27788559, 27788655, 27788724, 27788906, 27789265, 27789457, 27789813, 27789905, 27797438, 27798483, 27802101, 28467891, 28467923, 28467936, 28467976, 28567093, 28567168, 28873977, 28877281, 28877408, 28877602, 28877619, 28918563, 28918651, 28921455, 28951669, 28951709, 28951763, 28951804, 28951812, 29829333, 29829364, 29829883, 29829905, 29829996, 29830011, 29830746, 29831260, 29831467, 29835094.

## The eleven systems, and what is in each QM region

Read off `preprocess.results.qm_residues` in each directory's
`config.yaml` on 2026-11-01, not from the configs in this repo - the run
directories are the record of what was actually computed.

| directory | QM region | Arg90 | Lys226 |
| --- | --- | --- | --- |
| `CSR_crystal` | Zn, HOH485, His88, His164, Asp92, Arg90, Lys226, OH, mer | in | in |
| `CSR_crystal_v2` | same nine | in | in |
| `CSR_crystal_v3` | same nine | in | in |
| `CSR_gln93` | the nine **plus Gln93** | in | in |
| `CSR_hydrolysed` | Zn, mer, Asp92, His88, Arg90 - **no His164, no hydroxide, no HOH485** | in | **out** |
| `CSR_ligand` | Zn, OH, mer, Asp92, His88, Arg90, Lys226 - no His164, no HOH485 | in | in |
| `CSR_nac` | the nine | in | in |
| `CSR_ohH` | no preprocess results recorded | - | - |
| `CSR_unrestrained` | the nine | in | in |
| `CSR_unrestrained2` | the nine | in | in |
| `CSR_wt` | the nine | in | in |

**Arg90 and Lys226 are in the QM region by choice, not by detection.**
`preprocess.inputs.extra_qm_residues: {A: [90, 226], C: [1]}` puts them
there; neither is a zinc ligand, so the 2.5 A detection would never have
picked them up. That matters for the charge-deletion decomposition (¶153):
the two residues cannot be measured by deleting MM charges they do not
have, and **no equilibrated system we own leaves them both in MM.**
`CSR_hydrolysed` is the only one with Lys226 outside the QM region, and it
is the hydrolysed product with an incomplete first shell - the wrong state
and the wrong model. Measuring them that way needs a fresh `preprocess`
with a smaller `extra_qm_residues` and a re-equilibration of everything
downstream, which is a new model and not comparable with the numbers
already in the log.

## The structures that are live inputs today

| structure | what it is | made by |
| --- | --- | --- |
| `CSR_wt/protdft/scan_aspin/scan_00.pdb` | four-coordinate site, Asp92 out at 3.60 A, substrate present but not coordinated; the geometry both charge-deletion single points scored | 31143211 (¶140/¶143) |
| `CSR_wt/protdft/scan_aspin/scan_08.pdb` | the same system with Asp92 on the metal at 1.90 A, five ligands | 31143211 |
| `CSR_wt/protdft/scan_aspin/scan_08_A92_ASH.pdb` | scan_08 with Asp92 protonated (HD2 on OD2, anti, in plane) | `site.py protonate` (¶152) |
| `CSR_wt/protdft/qmmm_asph/QM_MM_optimized.pdb` | the four-coordinate sphere the neutral-Asp92 optimisation produced | 31186037 (¶153) |
| `CSR_wt/protdft/scan_michaelis2d/` | the 45-point 2D surface, all points converged | 31155441 (¶147/¶149) |
| `CSR_nac/md_naclong_results/` | the 500 ps NAC MD whose last frame fed `nacqm`, `nacarg` and `ohH` | 30421166 (¶101/¶102) |
| `CSR_wt/md_wt100b_results/` | 100 ns with the substrate tethered, first and last frames kept | 31016877 (¶138) |

## Trajectories deleted on 2026-11-03, and what is left in their place

Scratch was at 275/300 GB with `mich_eq5` queued and about 6 GB of trajectory
still to write, so ten `.dcd` files were deleted. Every one of them had its
`*_lastframe.pdb` (and `*_lastframe_imaged.pdb`) on disk first, which is what
was checked before anything was removed: the frames the decision log quotes are
all still here, only the frames between them are gone. 6 GB recovered,
269/300 GB after.

| deleted | why it was redundant |
| --- | --- |
| `CSR_wt/md_wt100c_results/NVTtrajectory.dcd` | the same trajectory as `NVTtrajectory_imaged.dcd`, which is kept and is what `postprocess` reads |
| `CSR_wt/md_wt100c_results/relaxbox_NPT{,_imaged}.dcd` | box equilibration; nothing is measured on it and its last frame is kept |
| `CSR_prod/md_prod100c_results/relaxbox_NPT{,_imaged}.dcd` | as above |
| `CSR_crystal/md_restrained/NVTtrajectory_imaged.dcd` | superseded by `CSR_crystal_v2` |
| `CSR_unrestrained/md_results/NVTtrajectory_imaged.dcd` | unrestrained, so invalid for coordination by the project's own rule |
| `CSR_unrestrained2/md_results/NVTtrajectory_imaged.dcd` | as above |
| `CSR_ligand/md_results/NVTtrajectory_imaged.dcd` | superseded build, Sep 24 |
| `CSR_hydrolysed/md_results/NVTtrajectory_imaged.dcd` | superseded build, Sep 24 |

Three were kept deliberately:

- `CSR_wt/md_wt100c_results/NVTtrajectory_imaged.dcd` - the 100 ns gate run.
- `CSR_prod/md_prod100c_results/NVTtrajectory.dcd` - the product arm's NaN at
  22.68 ns is still unexplained, and this run has **no** last frame precisely
  because it died, so the trajectory is the only record of it.
- `CSR_crystal_v2/md_results/NVTtrajectory_imaged.dcd` - the parent of
  `CSR_mich`, cheap at 0.29 GB.

So the ten ¶195 verdicts can be re-read but not re-measured frame by frame. If
a revision needs the distributions rather than the endpoints, those runs have
to be repeated, and that is a deliberate trade for the headroom.

### The second prune, later the same evening

`md_eq5` alone wrote **13 GB**, because `md` keeps every trajectory twice, raw
and imaged. Once ¶197's analysis had been collected, three more files went:

| deleted | size | why redundant |
| --- | --- | --- |
| `CSR/CSR_mich/md_eq5_results/NVTtrajectory.dcd` | 5.4 G | identical to `NVTtrajectory_imaged.dcd`, which is what `postprocess` actually read - `results.coordination_trajectory` names the imaged copy explicitly, so this was checked rather than assumed |
| `CSR/CSR_mich/md_eq5_results/relaxbox_NPT.dcd` | 550 M | box equilibration; nothing is measured on it |
| `CSR/CSR_mich/md_eq5_results/relaxbox_NPT_imaged.dcd` | 550 M | as above |

`md_eq5_results` went from 13 G to **5.7 G** and scratch from 286 to
**280/300 GB**. Kept: the imaged NVT trajectory, `NVTtrajectory_lastframe.pdb`
and `..._lastframe_imaged.pdb` - the frame every ¶197 number was measured on -
and `relaxbox_NPT_lastframe.pdb`. The three `warmup_MD_cycle*` pairs were left
alone, about 180 MB between them.

**The cost, stated because unlike the first prune this run backs a live open
question.** The 97.9% occupancy, the coordination distribution and the sphere
angles can all be recomputed from the imaged copy, so nothing in ¶197 is
unreproducible. What is gone is the unwrapped coordinates: any property that
imaging destroys - diffusion, or a distance measured across a periodic
boundary - would mean re-running the 10 ns.
