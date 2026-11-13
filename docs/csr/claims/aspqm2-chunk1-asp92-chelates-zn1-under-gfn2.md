---
type: observation
status: supported
evidence: ours
entries: [¶231, ¶232]
checked: [¶232]
---

# From a crystal-consistent start with the salt bridge intact, GFN2 brings Asp92 onto Zn1 as a bidentate chelate within 8 ps and keeps it there

`aspqm2` chunk 1, job 32416288, read live at 32 of 50 ps (32,122 xTB calls).
50 ps of GFN2-xTB/MM dynamics, **nothing restrained**, charge 2, QM region
`Zn, OH, mer, His88, His164, Asp92, Asn225, Arg90, Lys226`, started from
`md_relax10`'s last frame ([[relax10-holds-the-crystal-geometry-for-10-ns]]),
whose Asp92 sits where 6DQH puts it and whose Lys226–Asp92 salt bridge is
formed.

**The carboxylate crosses the coordination boundary between 7 and 8 ps and does
not come back for the remaining 24 ps.**

| frame (ps) | 1 | 4 | 7 | 8 | 9 | 15 | 20 | 25 | 32 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Zn–OD2 | 4.34 | 4.34 | 3.54 | 2.36 | 2.01 | 2.27 | 2.03 | 2.07 | 2.12 |
| Zn–OD1 | 5.90 | 5.59 | 4.74 | 3.52 | 3.41 | 2.31 | 2.28 | 2.23 | 2.29 |
| Zn–O2x | 2.50 | 2.12 | 2.21 | 2.51 | 2.60 | 4.83 | 4.94 | 6.31 | 5.81 |

Occupancy inside the 2.5 Å bond radius over all 32 frames: His88 NE2 **100%**
(2.09), His164 NE2 **100%** (2.12), hydroxide O01 **100%** (1.83), Asp92 OD2
**78%** (2.10), Asp92 OD1 **53%** (2.29), carbonyl O2x **19%** (2.16),
HOH3171 **3%**, HOH3095 **0%**. Coordination number 3 at 1 ps, 4 through
13 ps, **5 from 14 ps on**.

**It is a chelate, not a contact.** Zn–CG falls 5.50 → 2.52 Å (min 2.34), which
is the distance a symmetric bidentate carboxylate puts its own carbon at, and
both oxygens are inside the radius from 15 ps.

**Independently confirmed from a different file written by a different writer.**
`xtb_.xyz`, the live 130-atom QM region handed to xTB after frame 32, gives the
zinc's neighbours inside 2.9 Å as **O 1.813, N 2.023, O 2.030, O 2.074,
N 2.294, C 2.389** — two histidine nitrogens, the hydroxide, *two* carboxylate
oxygens and the carboxylate carbon. The DCD parse and the live QM geometry
agree without sharing any code.

**Three things did NOT happen, and each matters.**

- **Lys226 kept its proton.** All three NZ hydrogens sit at 0.98–1.09 Å for
  every frame, and O01 carries exactly one hydrogen (0.98–0.99 Å) throughout —
  it is still a hydroxide. So the collapse is **not** downstream of the Lys226
  proton transfer, which is how [[zinc-bound-crystal-holds-asp92-off-the-metal]]
  and ¶216 proposed to explain it away. See the caveat this puts on
  [[lys226-protonates-the-nucleophile-under-qmmm]].
- **The salt bridge held.** NZ–OD2 mean **2.87 Å** (min 2.53, max 4.24) against
  2.79 Å in 6DQH. The carboxylate is bonded to the metal at 2.10 Å *and*
  hydrogen-bonded to the ammonium at 2.87 Å at the same time, for 24 ps. Holding
  the lysine contact does not keep the carboxylate off the zinc.
- **The carboxylate arrived first and the substrate was displaced after.** Zn–OD2
  crosses inward between 7 and 8 ps; Zn–O2x crosses outward between 8 and 9 ps.
  So no vacancy opened first for the carboxylate to fill.

**What this run is not.** The Michaelis complex dissociates: with nothing
restrained the carbonyl leaves the metal (mean 4.32 Å) and the attack coordinate
O01–C4x goes from 3.4–4.5 Å through 22 ps to 5.6–6.2 Å after it. By criterion 6
this is not a Michaelis complex past ~22 ps, so **nothing energetic may be taken
from these frames** — the measurement here is about where the carboxylate goes,
which was the question asked.

## Supports
- [[gfn2-over-coordinates-zn1-from-a-crystal-consistent-start]]
