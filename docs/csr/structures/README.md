# Local copies of CSR coordinate files

`*.pdb` here is **gitignored** — these are multi-megabyte solvated frames that
belong on Bunya, where they are inputs to further jobs. This directory is for
having them to hand locally (PyMOL, measurement, figures), not for version
control.

## What is here, pulled through the code-server read route

| file | source on Bunya | what it is |
| --- | --- | --- |
| `l1_in_relax10.pdb` | written by `refstruct.py transfer`, ¶247 | 7UHI's moxalactam (`MX0`, 23 atoms) fitted onto `relax10`'s zinc site, RMSD 0.310 Å. The experimental Michaelis geometry **in our frame** — the reference ¶247 measured the 6.20 Å centroid separation against |

Regenerate `l1_in_relax10.pdb` with, from `CSR_mich/`:

    F=A:105:NE2=A:88:NE2,A:181:NE2=A:164:NE2,A:109:CG=A:92:CG
    python $R/refstruct.py transfer /scratch/user/uquadhik/refs/7UHI.pdb \
      md_relax10_results/NVTtrajectory_lastframe.pdb \
      --ligand=A:303:MX0 --zn=B:1:ZN --fit=$F --out=l1_in_relax10.pdb

The ligand residue name is **`MX0`**, not `MOX` — searching for the latter
finds nothing.

## How these got here: the remote, not a download

Both this machine and Bunya have the same `mcube` checkout and the same GitHub
remote, so a stripped frame travels by `git push` on Bunya and `git pull` here.
That is the whole mechanism. `*_nowat.pdb` is **tracked** for exactly this
reason; full solvated frames, trajectories and `.csv` are not.

`.gitattributes` marks `*.pdb`, `*.xyz` and `*.cif` as `-text` so no
line-ending conversion happens and a round-tripped frame is byte-identical.

Each file was checked after arriving by recomputing distances the decision log
already records — `relax10` Zn-OD2 5.153 / Zn-OD1 4.327 / Zn-CG 5.192 /
Zn-NZ 5.276, `aspxtb8` 5.092 / 5.830 / 5.777 / 5.044, `prodpull` 4.465 / 4.852
/ 5.131 / 5.474 — all reproduced to the third decimal. That is also the test
that `strip_water.py` is faithful.

## The full solvated frames, and how to get them

There is **no SSH private key on this machine** — `~/.ssh/` holds only
`known_hosts`, with `bunya.rcc.uq.edu.au` in it. So `scp` cannot authenticate
non-interactively, and the code-server read route returns file contents through
the assistant's context, which a 4 MB frame will not fit in usefully.

Run `pull-from-bunya.sh` yourself to fetch them; it will ask for your Bunya
credentials, which is the step that cannot be automated from here.
