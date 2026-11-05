---
type: observation
status: supported
evidence: ours
entries: [D1, D17, D144, D207]
checked: [D144]
---

# No reaction coordinate has ever been computed in this project

`neb` has never completed, `ts` has never run: there is no barrier, no
transition state, no IRC, no intermediate (audit headline, confirmed against
D144). CI-NEB was never deselected - it is in D1, the first entry, and survived
every revision of the plan - and the machinery is working. `neb_step1` names
`protdft/scan_attack/scan_15.pdb` as its product and **that file has never
existed**, because the scan was never run (D207).

**The diagnosis worth keeping** (D207): each blocker was investigated on its
own terms rather than as a question about the endpoint. Asp92's coordination
got three campaigns, the protonation state got three, and the loop could not
terminate because every attempt to settle the reactant used our own minima as
the evidence - which D206 rules inadmissible for exactly that. **A question
whose only admissible evidence class cannot settle it will reopen
indefinitely**, and it did.

**And a correction to how the number should be framed.** NEB was our
substitution for Imhof's path network, not the field's method. No MBL mechanism
study in this project's survey uses NEB; CphA's ~14 kcal/mol - the mono-zinc
number we would compare against - is a **PMF** barrier, and D32 was a user
decision to replicate that. So a minimum-energy barrier quoted beside it
compares two different quantities (D207, as qualified by the user the same day).

## Supports
- [[thesis-mechanism-differs-from-canonical-b3]]

## Refutes
- [[rqk-cations-replace-zn2-electrostatically]]
