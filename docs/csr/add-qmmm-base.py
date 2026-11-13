# -*- coding: utf-8 -*-
"""Give a run directory the base `qmmm:` section its named sections inherit.

`qmmm` layers `qmmm.inputs` under `<section>.inputs`, exactly as `md` does, and
refuses to run without the base: `KeyError: No qmmm.inputs section in
config.yaml`. A directory built by `build-mich.py` has only `preprocess:` and
`md:`, because that is all the two stages it ran needed - so the first QM/MM
section submitted against it dies in 5 seconds (job 31738709, ¶198).

This writes that base. It is a script rather than a `site.yaml` section because
`sections.py` writes `<stage>_<name>` and **refuses a name equal to its own
stage**, so a bare `qmmm:` key cannot be produced that way.

**The MM side is taken from the run's own `md` inputs, not from a template.**
The force fields a QM/MM run uses must be the ones its MD used, or the MM region
is a different system from the one that was equilibrated and nothing downstream
would notice: the energies are all relative, the geometry looks normal, and the
only symptom is a number that cannot be compared with anything. So
`default_forcefields` comes from `md.inputs` and the ligand XML from
`preprocess.results.ligand_xml`, and every file in the list is checked to exist
before anything is written.

The QM/MM-only values follow `CSR_wt`'s base, which every QM/MM run in the
project has used: 1 fs timestep, 300 K, an 8 A active region centred on the
metal, `opt_maxiter` 300. `charge` is deliberately left at the value a caller
must override - each section states its own (¶152), and a base charge that
happens to suit one arm is how an SCF converges on the wrong species while
looking entirely normal.

Usage, from the run directory:

    python add-qmmm-base.py                 # config.yaml here
    python add-qmmm-base.py path/to/config.yaml
    python add-qmmm-base.py --force         # overwrite an existing qmmm.inputs
"""

import sys
from pathlib import Path

import yaml

# Values that are properties of the stage rather than of a particular arm.
QMMM_DEFAULTS = {
    "cores": 4,
    "temp": 300,
    "timestep": 0.001,
    "time": 50,
    "mult": 1,
    "implicit": False,
    "skip_equilibrate": True,
    "skip_optimize": False,
    "singlepoint": False,
    "actradius": 8,
    "opt_maxiter": 300,
}

# The active region is centred on the metal, so this one IS system-specific.
# Every CSR run uses B:1:ZN; another system must pass `--origin`.
DEFAULT_ORIGIN = "B:1:ZN"


def main(argv):
    force = "--force" in argv
    args = [a for a in argv[1:] if not a.startswith("--")]
    config = Path(args[0]) if args else Path("config.yaml")
    if not config.is_file():
        raise SystemExit(f"no such config: {config.resolve()}")

    text = config.read_text(encoding="utf-8")
    doc = yaml.safe_load(text) or {}

    if (doc.get("qmmm") or {}).get("inputs") and not force:
        existing = doc["qmmm"]["inputs"]
        raise SystemExit(
            f"{config} already has a qmmm.inputs with {len(existing)} keys, and "
            "overwriting it would change the MM treatment of every QM/MM "
            "section in this run. Pass --force only if that is what you mean.")

    md_inputs = (doc.get("md") or {}).get("inputs") or {}
    pre_results = (doc.get("preprocess") or {}).get("results") or {}

    forcefields = list(md_inputs.get("default_forcefields") or [])
    if not forcefields:
        raise SystemExit(
            f"{config} has no md.inputs.default_forcefields, so there is "
            "nothing to copy the MM treatment from. Writing a guessed list "
            "would give the QM/MM a different MM system from the MD's.")

    ligand_xml = pre_results.get("ligand_xml") or md_inputs.get("ligand_xml")
    if ligand_xml:
        # Keep it relative to the run directory where it already is, so the
        # entry survives the directory being moved (the CSR/ move of 2026-11-03
        # rewrote absolute paths once already).
        here = config.resolve().parent
        try:
            ligand_xml = Path(ligand_xml).resolve().relative_to(here).as_posix()
        except ValueError:
            ligand_xml = Path(ligand_xml).as_posix()
        forcefields.append(ligand_xml)
    else:
        print("[WARNING] no preprocess.results.ligand_xml - the force-field "
              "list will have no ligand parameters. Correct only if this "
              "system has no non-standard residue.")

    missing = [f for f in forcefields
               if "/" in f and not f.startswith("amber")
               and not (config.resolve().parent / f).is_file()]
    if missing:
        raise SystemExit(
            "these force-field files are named but not present: " +
            ", ".join(missing) +
            f"\n(resolved against {config.resolve().parent})\nNothing written.")

    origin = next((a.split("=", 1)[1] for a in argv if a.startswith("--origin=")),
                  DEFAULT_ORIGIN)

    inputs = {"forcefields": forcefields}
    inputs.update(QMMM_DEFAULTS)
    inputs["actregion_origin"] = origin

    # Read-modify-write under the project's lock: a stage writing its results
    # back into this same file at the same moment would otherwise lose one of
    # the two edits, and the loser is silent.
    def edit(current):
        current.setdefault("qmmm", {})["inputs"] = inputs
        return current

    try:
        from mcube.preprocess import _locked
    except Exception:
        config.write_text(
            yaml.safe_dump(edit(doc), sort_keys=False, default_flow_style=False),
            encoding="utf-8")
        print("[WARNING] mcube not importable, so this was written WITHOUT the "
              "config lock. Safe only because no job is writing this config.")
    else:
        # `_locked` takes a Path and derives the lock file with `with_name`, so
        # a str gets an AttributeError three frames down.
        with _locked(config.resolve()):
            fresh = yaml.safe_load(config.read_text(encoding="utf-8")) or {}
            config.write_text(
                yaml.safe_dump(edit(fresh), sort_keys=False,
                               default_flow_style=False),
                encoding="utf-8")

    print(f"[OK] wrote qmmm.inputs into {config}")
    for key in ("forcefields", "timestep", "time", "actradius", "actregion_origin"):
        print(f"[INFO]   {key}: {inputs[key]}")
    print("[INFO] `charge` is NOT set here on purpose: every section states its "
          "own (¶152).")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
