---
type: claim
status: open
entries: [D197, D199, D200, D211]
checked: [D197, D199]
---

# Asp92 coordinating Zn1 interferes with hydrolysis

The mechanistic half of the thesis: it is what makes restraining the carboxylate
worth doing, rather than a structural curiosity.

**The argument.** Zn1 must present a hydroxide to the substrate carbonyl at an
attacking angle. A five-coordinate Zn1 cannot
([[a-fifth-zinc-donor-closes-the-attack-angle]]): across three frames the attack
angle ran 39/56/68 deg against a 95-120 deg requirement, and the one frame where
the offending ligand left was rescued by a different fifth donor within the same
frame. Asp92 on the metal is such a fifth donor. In the QM/MM lineage where it is
on the metal, **Zn-Asp92 OD2 averages 1.93 A over 50 frames, 100% inside 2.5 A,
and the nucleophile is never poised** (O01-C4x mean 3.81 A, D211).

**Why this is `open` and not `supported`, in three parts.**

1. **The crowding was never measured with Asp92 as the fifth donor.** It was
   measured with the substrate's O1x and with Asn225. The inference that a
   carboxylate oxygen behaves the same way is reasonable and untested.
2. **The one frame that made the crowding look unavoidable has been
   overturned.** D198 argued the site keeps a fifth donor whatever is done,
   because Asn225 OD1 took the vacated site at 2.19 A; **D200 measured Asn225
   leaving**, to 4.18 A closest and 5.86 A mean. So a four-coordinate sphere is
   not ruled out by anything measured, and the crowding may be an MM property
   rather than a property of the site.
3. **All of it is `ours`, and most of it is MM.** D211's own finding was that
   O1x-on-metal is an MM artefact that QM/MM does not reproduce. The same
   suspicion applies in reverse here: the Asp92 collapse is seen in QM/MM but in
   a lineage that started from an already-collapsed frame and had Lys226
   deprotonated, so one event had already happened at t=0.
4. **The literature reads the aspartate differently.**
   [[asp120-mutants-lose-activity-in-b1-and-b3]]: Asp120 mutants lose activity in
   every subclass tested, and the published explanation is that the residue
   positions **Zn2**, not that it must be kept off Zn1. For an enzyme with no Zn2
   that explanation does not transfer, but neither does it support us.

**What would settle it:** a barrier computed with Asp92 restrained off the metal
and again with it on, in the same system at the same level.
[[question-does-restraining-asp92-lower-the-barrier]].

## Requires
- [[question-does-restraining-asp92-lower-the-barrier]]

## Supports
- [[thesis-lys226-restrains-asp92-in-place-of-zn2]]

## Supported by
- [[a-fifth-zinc-donor-closes-the-attack-angle]]
- [[our-models-put-asp92-on-the-zinc]]
