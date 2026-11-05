# Working in an argument map - operating instructions

For an agent or a person arriving at a directory of argmap notes. Read
`SPEC.md` for the format; this file is what to *do*.

---

## Before you touch anything

1. **Run the validator.** `python argmap.py <dir>` (or `python check.py` if the
   project ships a shim). Its findings are the current state of the argument's
   health, and you should not add to a map whose existing findings you have not
   read.
2. **Read the thesis node.** It tells you what the map is for. If what you are
   about to add does not bear on it, it does not belong here.
3. **Find the project's maintained-state file and its retractions index** if it
   has them (`argmap.json` names them). A journal's headline is not a current
   value.

## Starting a map from scratch

1. **Get the thesis as an assertion.** If the project's framing is a question
   ("what is X's mechanism?"), it cannot be mapped - nothing supports or refutes
   a question. Ask for the claim the work is making. This is the single most
   important step and the easiest to skip.
2. **Write the thesis node first**, with its `requires` edges: the premises
   that must hold for the thesis to be available at all. These are usually
   fewer than expected, and listing them is where unsupported premises surface.
3. **Work downward from the thesis, not forward from the journal.** Reading a
   journal front to back produces a chronology with node headings. Start from
   "what would have to be true" and go looking for each piece.
4. **Target tens of nodes, not hundreds.** See SPEC's Scope section.

## Adding a node

1. **Decide the type.** If it is a datum, it is an `observation` and needs an
   evidence class. If it is something the project asserts on the strength of
   data, it is a `claim`. If it asserts nothing, it is a `question`.
2. **Name the file for the claim, not the topic.** `asp92-is-not-a-zn1-ligand`,
   not `asp92`. A node whose name is a noun will accumulate unrelated
   assertions and cannot be refuted.
3. **State it in the heading, in one line.** The heading is the claim.
4. **In the body, give the numbers and cite where each came from.** A claim
   that cannot be traced to a measurement, a source, or a named calculation is
   an assertion - mark it `status: unsupported` and say so in the body. Do not
   promote it by writing confidently.
5. **Add the edges on this node** (see below), then run the validator.

## Adding an edge

Write it on the node that **makes** the assertion. The reverse index is built
for you, so you never need to edit a second file to make an edge visible - but
do add a mirror heading where a reader of the other node needs to see it, and
always on a refuted claim.

- `supports` / `refutes`: between a node and the thing it is evidence for or
  against. Only add `supports` from an observation when the observation really
  is evidence for that claim, not merely adjacent to it.
- `requires`: for a premise, and for an open `question` that gates a claim.
- `superseded-by`: see below.

**If you find yourself wanting a fifth edge type, you want a sentence in the
body instead.** Nuance goes in prose; the edges are what gets checked.

## The discipline rules

**Never silently delete or quietly rewrite a node that turns out to be wrong.**
Set its `status` to `void`, add `## Superseded by` pointing at what replaced
it, and leave the body saying what was believed and why. A map that erases its
mistakes cannot show you that you are making one again - which is the whole
reason the map exists beside an append-only journal rather than instead of one.

**Do not mark something `supported` to make a check pass.** The checks exist to
find exactly that move. If check 2 flags a claim, the options are: demote the
claim, refute the refuter with evidence, or record in the body why the
refutation does not bite. Flipping the status is not one of them.

**Do not resolve a finding by deleting an edge.** A refutation you remove is a
contradiction you have hidden.

**An `ours` observation cannot promote a geometry, a structure or a fact it was
handed.** A calculation performed on an input is evidence about the
calculation. This is why the evidence classes exist, and check 1 is its
enforcement.

**Cite the number, not the conclusion.** "4.20 A in 6NC5 (D208)" survives a
retraction; "Asp92 is far from the zinc" does not travel with its provenance
and will be quoted after the measurement is superseded.

**When you cite an entry the retractions index lists as overturned**, read the
superseding entry first, make sure the body reflects it, then add the entry id
to the node's `checked` list. Leave it unacknowledged only when the node's
*point* is that its own support is retracted.

**Edit the project's maintained-state file in the same change**, if it has one
and the node changes a settled value. The graph is not a substitute for it.

## Finishing

Run the validator again and **read every finding**. A finding you cannot
resolve is a result: record it where the project records results. Do not leave
the map in a state whose findings you have not either fixed or explained.
