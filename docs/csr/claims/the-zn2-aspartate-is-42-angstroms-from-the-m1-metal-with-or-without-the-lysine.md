---
type: observation
status: supported
evidence: experimental
entries: [¶285]
checked: [¶284]
---

# The Zn2 aspartate sits ~4.2 angstroms from the M1 metal in every structure that has one, whether or not the lysine is present

Measured **directly, with no superposition**, in every available structure that
actually contains metals. Numbering is each file's own.

| structure | lysine? | metal | Asp->metal (closest OD) |
| --- | --- | --- | --- |
| 6NC5 ZN302 (M1) | **no** (225 is HIS) | 2 His + 2 HOH | Asp91 OD2 **4.20** |
| 6NC5 ZN301 (M2) | no | Asp,His,His + HOH | Asp91 OD2 **1.95** |
| `refine_16` ZN B2 (M1) | **no** (212 is HIS) | 2 His + 3 HOH | Asp78 OD1 **4.12** |
| `refine_16` ZN B1 (M2) | no | Asp,His,His + 3 HOH | Asp78 OD2 **2.06** |
| `refine_12` ZN B1 (M2) | no | Asp,His,His + 3 HOH | Asp78 OD2 **2.13** |
| 6DQH ZN A401 (M1) | **yes** (Lys260) | **3 His** + PO4 + 2 HOH | Asp120 OD1 **4.28** |

**4.20 and 4.12 without the lysine; 4.28 with it.** The aspartate is never an M1
ligand, and the lysine makes no measurable difference to how far out it sits.
Where it does coordinate a metal it coordinates **its own**, at the M2 site,
1.95-2.13 Å - it recruits a second metal rather than reaching across for the
first.

**This retracts a number.** ¶284 reported the aspartate at **3.00 Å** from the
M1 metal in wild type, obtained by superposing metal-bearing 6NC5 onto
metal-free 6DN4 and carrying a metal across. 6NC5 **has its own metals**, so the
fit was never needed, and measured in place the distance is **4.20 Å**. The
3.00 Å figure was an artefact of placing a metal into an apo structure - the
same construction ¶184 was voided for - and is withdrawn, together with the
conclusion ¶284 drew from it, that the aspartate is *closer* when the lysine is
present.

**¶284 also called `refine_16` "the one structure in the world with a
CSR-1-family M1 metal and no lysine".** It is not. **6NC5** has an M1 metal
(ZN302: His87, His163 and two waters), position 225 is a histidine, and HOH493
bridges it to the M2 metal. That is a second, independent, lysine-free
observation, and it agrees with ours.

## Supports
- [[question-how-far-is-asp92-from-zn1-in-wild-type-csr1]]
- [[asp92-is-never-a-zn1-ligand-in-b3-structures]]
