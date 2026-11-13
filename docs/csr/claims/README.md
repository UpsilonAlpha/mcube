# CSR-1: the argument map

The claims that bear on **one thesis**: *Lys226 replaces Zn2 in the CSR-1
wild-type mechanism by restraining the bridging aspartate, allowing MBL
activity* ([[thesis-lys226-restrains-asp92-in-place-of-zn2]]).

**The thesis was sharpened on 2026-11-06** from "CSR-1's wild-type mechanism is
different from the canonical B3 mechanism", which is kept as a claim
([[thesis-mechanism-differs-from-canonical-b3]]) because the new thesis entails
it. The old one is a negative claim and cannot be the content of a paper on its
own; the positive claim about the lysine can, and it names a mechanism, so its
premises can be listed and checked - which is what the map is for.

One markdown file per node, typed links between them, and a validator that
checks the result. The format is generic: `../../argmap/GUIDE.md` describes the
whole system in one file, `../../argmap/SPEC.md` is the terse format reference,
and `../../argmap/AGENT.md` is how to work in it. This directory is the
instance.

## The one rule

**The graph holds what we believe; the log holds what happened.**

`../../csr-decision-log.md` is not migrated and is never edited to match this
directory. Its chronology is evidence: ¶178 timestamps a blind prediction
precisely so it is checkable later, and "this contradiction sat unnoticed for
eight days" is a fact about the process. Every node here cites the log entries
it came from; for a *current* value, `../state.md` remains the authority.

**What this layer adds is linkage.** The log's failure mode is not search: ¶58
asserted a conclusion, ¶59 retracted it and then re-asserted it in its own
closing line, and the retracted half was quoted as the project's position for
two weeks. `../retractions.md` patches that for a reader who checks the index.
An inbound `refutes` edge patches it for a reader who does not.

The live example, which the validator finds in one command: ¶201 closed the
CphA route for the third time because "Asp92 is a zinc ligand in every state
this model reaches", and that premise is refuted by measurement
([[asp92-is-never-a-zn1-ligand-in-b3-structures]]). See
[[cpha-route-is-closed]] - **demoted to `open` on 2026-11-06** as a result of
exactly that edge.

Three further nodes moved the same day, all found by running the validator
rather than by reading: [[asp92-is-a-zn1-ligand-in-our-models]] was `supported`
while its own body said it is "wrong by 2 A or more" against the crystals (now
`void`); [[mono-zinc-mbl-precedents-are-zn2-enzymes]] was `supported` on a
GOB-18 premise its own body says the deposit does not evidence (now `open`); and
[[a-fifth-zinc-donor-closes-the-attack-angle]] had to be rewritten when the
index showed ¶200 overturning the Asn225 frame that made its generalisation
look safe.

## Four edge types

`supports`, `refutes`, `requires`, `superseded-by`. Written as body sections on
the node that makes the assertion; the validator builds the reverse index, so
`## Supports` and `## Supported by` produce the same edge. **Four, and a fifth
is not to be added** - link discipline is what this project has failed at.

## Three evidence classes

Every `observation` carries one: **`experimental`** (a deposited structure, a
measurement from one) > **`literature`** (a published study) > **`ours`** (our
own calculation). This is ¶206's hierarchy as data instead of a rule to
remember, and it has a mechanical consequence: **a claim whose only support is
`ours` is weakly held by construction**, and the validator says so.

## Running the validator

```
python check.py
```

From this directory, no arguments, no dependencies. It prints nine numbered
sections: note format (0) and one-thesis (0b); what each supported claim rests
on, following `supports` transitively (1); supported claims with a standing
refutation (2); void nodes with no replacement and live nodes still resting on
one (2b); premises the thesis requires, transitively, with no observation
beneath them (3); supported nodes citing entries `../retractions.md` lists as
overturned (4); broken wikilinks (5a); and nodes with no path to the thesis
(5b). Every finding carries a `fix:` line. Exit 1 if anything is flagged, 2 on
a structural error.

`python check.py --json` gives the same report as JSON for a tool or an agent
to act on, and `--strict` counts the advisory notes towards the exit status.

`argmap.json` points it at this project's log and retractions index.

## Adding a node

1. Run `python check.py` and read the existing findings first.
2. One file, kebab-case, named for the **claim** and not the topic.
3. Frontmatter: `type` (`thesis`/`claim`/`observation`/`question`), `status`
   (`supported`/`open`/`unsupported`/`withheld`/`void`), `evidence` on
   observations only, `entries` for the log citations - required on a
   `supported` node, and all five keys are validated, so a typo is a finding
   rather than a silent exemption from scrutiny.
4. Body: the claim as the heading, then the numbers and the entry each came
   from. **A claim that cannot be traced to a measurement, a paper, or a named
   calculation is an assertion - mark it `status: unsupported`.**
5. Edges as `## Supports` / `## Refutes` / `## Requires` / `## Superseded by`.
6. Re-run the validator. Never resolve a finding by flipping a status or
   deleting an edge; a node that turns out wrong goes to `status: void` with a
   `superseded-by` edge, and stays.

## Scope

~38 nodes, not 208. Operational history - job ids, Slurm rates, scratch quota,
config bugs, resume wiring - stays in the log: it has no argumentative role and
converting it would only add noise.

Deliberately withheld and visible as such:
[[mutant-and-double-mutant-kinetics]]. **Do not go looking for the kinetics**
(¶178).
