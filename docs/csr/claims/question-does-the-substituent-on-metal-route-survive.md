---
type: question
status: open
entries: [¶126, ¶196, ¶278]
---

# Does the substituent-on-metal arrangement survive as a mechanism, or was it an artefact of the docked pose?

**This hypothesis has been in `site.yaml` since ¶126 with its own full criteria
set, `substituent_on_metal`, and has never had a node on this map.** Converting
the gates found it: a mechanism the pipeline can accept a structure under, with
nothing here to argue about it.

**What it is.** The grafted pose-2 structure presents the metal holding the
6-alpha-hydroxyethyl **O1x at 2.66 Å** while the beta-lactam carbonyl points
away at **4.50 Å**, and the zinc-bound hydroxide sits on the attack line at
**3.17 Å and 106 degrees** - inside the Burgi-Dunitz window ¶34 never reached
in 791 frames - with Asp92 OD1 at 2.87 Å.

It is a **hybrid** of the two routes ¶124 separated: the nucleophile is
metal-bound and on the attack line, as in thermolysin and carboxypeptidase A,
while the metal holds a substituent oxygen rather than the carbonyl, as GOB-18
does with the C-3 carboxylate.

**Why it stays a question rather than becoming a claim.** Nothing says a
6-alpha-hydroxyethyl oxygen belongs on a Zn1-site metal, and `site.yaml` marks
its own window honestly - `Zn - O1x, range [2.0, 2.8], source: NOT a reference
value - ¶126's own 2.66 Å plus tolerance`. So the gate for this hypothesis is
calibrated on the single structure that suggested it, which is circular unless
an independent structure or a barrier separates the two routes.

**It also has a live claim on the attack geometry that the preferred route does
not.** ¶126's 106 degrees was measured with **nothing restraining that
coordinate**, where the carbonyl-on-metal arm has never reached the window
without one. That is a real argument and the reason this is not simply dropped.

**What would settle it:** a barrier from each arrangement on the same topology,
or an experimental complex of a carbapenem with a 6-alpha-hydroxyethyl group on
any B3 enzyme showing which oxygen the metal takes.

## Supports
- [[every-acceptance-gate-encodes-a-mechanism]]

## Refutes
- [[carbonyl-coordinates-zn1-in-the-michaelis-complex]]
