---
type: observation
status: supported
evidence: literature
entries: [¶277]
checked: [¶277]
---

# S214D changes two things at once, and the clean S-to-A test says the serine is not required

**The substitution is confounded by construction.** It removes the intersite
serine *and* installs a carboxylate at the same position, adding sterics and a
negative charge a few Å from the nucleophile. No result from it can separate
"the serine was doing something" from "the aspartate is doing something".

**The manuscript's own argument runs against the serine being required.** It
cites Hou et al. 2017, who made the clean **S-to-A** substitution in AIM-1: the
enzyme "was still active against penicillins, albeit with reduced turnover and
affinity (kcat = 251 s-1, KM = 434 uM), thus suggesting the serine residue is
not strictly required". S-to-E was far worse (kcat 0.4, KM 57,000 uM),
attributed to "the steric interference of the bulkier residue". The paper reads
the aspartate result as "likely involves steric hindrance, but ... activity was
completely abolished suggests further interactions".

**And the phenotype is not one thing**, so there is no single observable to
reproduce: penicillins abolished ("no activity was detected in aspartate mutants
against any of the penicillin substrates"); carbapenems - our substrate - down
"3-5 orders of magnitude" in turnover with "only ... mild changes to binding
affinity"; cefuroxime abolished; cefaclor and cephalothin showing "no
substantial difference in kcat or KM". The abstract's "abolishes the
antibiotic-degrading function" holds for penicillins and cefuroxime only.

So the only computational observable S214D offers is a **barrier** - the most
expensive quantity in this project and the one most exposed to the unresolved
starting-frame question ([[question-is-the-reactant-pose-experiment-anchored]]).
It keeps a cheaper **structural** arm for later: whether the aspartate
carboxylate takes the hydroxide's place or chelates the metal, which a geometry
optimisation answers.

Quoted for interpretation only. The kinetics do not enter any barrier or
calibrate any simulation.

## Supports
- [[the-plain-double-mutant-is-the-mutant-arm]]
