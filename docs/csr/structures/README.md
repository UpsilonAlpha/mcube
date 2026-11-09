# Local copies of CSR coordinate files

`*.pdb` here is **gitignored** — these are multi-megabyte solvated frames that
belong on Bunya, where they are inputs to further jobs. This directory is for
having them to hand locally (PyMOL, measurement, figures), not for version
control.

## What is here, pulled through the code-server read route

| file | source on Bunya | what it is |
| --- | --- | --- |
| `l1_in_relax10.pdb` | written by `refstruct.py transfer`, D247 | 7UHI's moxalactam (`MX0`, 23 atoms) fitted onto `relax10`'s zinc site, RMSD 0.310 Å. The experimental Michaelis geometry **in our frame** — the reference D247 measured the 6.20 Å centroid separation against |

Regenerate `l1_in_relax10.pdb` with, from `CSR_mich/`:

    F=A:105:NE2=A:88:NE2,A:181:NE2=A:164:NE2,A:109:CG=A:92:CG
    python $R/refstruct.py transfer /scratch/user/uquadhik/refs/7UHI.pdb \
      md_relax10_results/NVTtrajectory_lastframe.pdb \
      --ligand=A:303:MX0 --zn=B:1:ZN --fit=$F --out=l1_in_relax10.pdb

The ligand residue name is **`MX0`**, not `MOX` — searching for the latter
finds nothing.

## Why the large files are not here, and how to get them

There is **no SSH private key on this machine** — `~/.ssh/` holds only
`known_hosts`, with `bunya.rcc.uq.edu.au` in it. So `scp` cannot authenticate
non-interactively, and the code-server read route returns file contents through
the assistant's context, which a 4 MB frame will not fit in usefully.

Run `pull-from-bunya.sh` yourself to fetch them; it will ask for your Bunya
credentials, which is the step that cannot be automated from here.
