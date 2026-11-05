# CSR-1: the argument map

The claims that bear on **one thesis**: *CSR-1's wild-type mechanism is
different from the canonical B3 metallo-beta-lactamase mechanism*
([[thesis-mechanism-differs-from-canonical-b3]]).

One markdown file per node, typed links between them, and a validator that
checks the result. The format is generic and documented in
`../../argmap/SPEC.md`; how to work in it is `../../argmap/AGENT.md`. This
directory is the instance.

## The one rule

**The graph holds what we believe; the log holds what happened.**

`../../csr-decision-log.md` is not migrated and is never edited to match this
directory. Its chronology is evidence: D178 timestamps a blind prediction
precisely so it is checkable later, and "this contradiction sat unnoticed for
eight days" is a fact about the process. Every node here cites the log entries
it came from; for a *current* value, `../state.md` remains the authority.

**What this layer adds is linkage.** The log's failure mode is not search: D58
asserted a conclusion, D59 retracted it and then re-asserted it in its own
closing line, and the retracted half was quoted as the project's position for
two weeks. `../retractions.md` patches that for a reader who checks the index.
An inbound `refutes` edge patches it for a reader who does not.

The live example, which the validator finds in one command: D201 closed the
CphA route for the third time because "Asp92 is a zinc ligand in every state
this model reaches", and that premise is refuted by measurement
([[asp92-is-never-a-zn1-ligand-in-b3-structures]]). See
[[cpha-route-is-closed]].

## Four edge types

`supports`, `refutes`, `requires`, `superseded-by`. Written as body sections on
the node that makes the assertion; the validator builds the reverse index, so
`## Supports` and `## Supported by` produce the same edge. **Four, and a fifth
is not to be added** - link discipline is what this project has failed at.

## Three evidence classes

Every `observation` carries one: **`experimental`** (a deposited structure, a
measurement from one) > **`literature`** (a published study) > **`ours`** (our
own calculation). This is D206's hierarchy as data instead of a rule to
remember, and it has a mechanical consequence: **a claim whose only support is
`ours` is weakly held by construction**, and the validator says so.

## Running the validator

```
python check.py
```

From this directory, no arguments, no dependencies. It reports: claims held on
our own calculations alone; supported claims with a standing refutation;
premises the thesis requires with no support; supported nodes citing entries
`../retractions.md` lists as overturned; broken wikilinks; and nodes with no
path to the thesis. Exit 1 if anything is flagged.

`argmap.json` points it at this project's log and retractions index.

## Adding a node

1. Run `python check.py` and read the existing findings first.
2. One file, kebab-case, named for the **claim** and not the topic.
3. Frontmatter: `type` (`thesis`/`claim`/`observation`/`question`), `status`
   (`supported`/`open`/`unsupported`/`withheld`/`void`), `evidence` on
   observations only, `entries` for the log citations.
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
(D178).
