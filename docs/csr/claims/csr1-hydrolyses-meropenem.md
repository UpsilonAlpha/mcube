---
type: claim
status: open
entries: [D22, D40, D178]
---

# CSR-1 hydrolyses meropenem

The substrate premise of the whole project, and **nothing in the record
measures it.** D22 is explicit: there are **no published kinetics for CSR-1
with meropenem** and no published mechanistic study. Pedroso 2020 characterises
CSR-1 as a broad-spectrum antibiotic-degrading B3-RQK MBL, which is not a
turnover number for this drug.

The measurement exists and is deliberately withheld
([[mutant-and-double-mutant-kinetics]]) under the blind-prediction protocol
(D178). So this node is `open` by design rather than by neglect, and its
support set in the graph is **empty on purpose**: adding the withheld kinetics
as support here would be claiming evidence nobody in this project has read.

**Why it is worth a node.** A referee will ask it, and the honest answer is
that the computed mechanism is being built for a reaction whose rate the
computation has not been shown. Flag it; do not quietly assume it.

## Supports
- [[thesis-mechanism-differs-from-canonical-b3]]
