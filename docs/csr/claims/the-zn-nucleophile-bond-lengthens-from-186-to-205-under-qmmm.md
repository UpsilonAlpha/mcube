---
type: observation
status: supported
evidence: ours
entries: [¶48, ¶118]
checked: [¶48, ¶118, ¶278]
---

# The Zn-hydroxide bond measures 1.83-1.86 angstroms under MM and relaxes to 1.99-2.05 under QM/MM

Measured twice independently, in ¶48 and again in ¶118: the same coordinate
reads **1.83-1.86 Å** on an MM-equilibrated frame and **1.99-2.05 Å** once
QM/MM takes over, with no restraint on it in either case.

The rest of the sphere is off in the other direction at the same time -
hydroxide short at 1.86, **His88 2.22** and **His164 2.32** long - so the error
is not a uniform scale factor but a systematic distortion of the Li/Merz 12-6
site. Crystal values for comparison: thermolysin 2.01, carboxypeptidase A 2.19.

So an absolute Zn-donor bond length measured on an MM frame reports the force
field. The matched DFT-versus-GFN2 pair puts the same bond at **MM 2.035 /
GFN2 1.776 / DFT 1.910** - GFN2 0.134 Å short of DFT, histidines 0.073-0.097
short ([[gfn2-overstates-the-metal-coordination-preference]]), so even between
two quantum levels the number moves by more than the width of the gate.

## Supports
- [[mm-cannot-be-gated-on-absolute-zn-donor-bond-lengths]]
