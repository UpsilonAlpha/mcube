"""docs/csr/site.py protonate mode: pytest tests/test_protonate.py

Asp92 is modelled deprotonated and coordinating the zinc. Testing the
alternative - a neutral aspartate, which would not coordinate - is blocked on
nothing but structure: the PDB carries no proton on the carboxylate and the
residue is named ASP, so ff14SB's ASH template does not match and the run dies
in `createSystem` rather than saying what is wrong.

What is worth pinning is the geometry and the refusals. A proton 0.98 A from the
oxygen but out of the carboxylate plane, or placed syn instead of anti, is a
structure that looks fine in a viewer and spends the first picoseconds of QM/MM
undoing a contact - and the comparison it is built for rests on the two states
differing by one hydrogen and nothing else.

The synthetic structures here are built in the test rather than read from a real
one, so the carboxylate geometry is known exactly and the assertions are about
the construction rather than about one frame. Needs pyyaml; the force-field
check needs OpenMM (the `mcube-test` env) and is skipped without it. ASH itself
is not installed anywhere and is not needed.
"""

import importlib.util
import math
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "docs" / "csr" / "site.py"

# The carboxylate: C-O 1.25 A, O-C-O 123 degrees, built in the z = 0.8 plane so
# "in the plane of CG/OD1/OD2" is a real constraint and not satisfied by accident.
CG = (0.200, 1.200, 0.800)
HALF_OCO = math.radians(123.0 / 2.0)
OD1 = (CG[0] + 1.25 * math.cos(HALF_OCO), CG[1] + 1.25 * math.sin(HALF_OCO), CG[2])
OD2 = (CG[0] + 1.25 * math.cos(-HALF_OCO), CG[1] + 1.25 * math.sin(-HALF_OCO), CG[2])


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("site_protonate", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def run(mod):
    """`protonate` with --to/--out passed the way the command line passes them.

    OPTS is module state, so it is cleared afterwards: a leaked --out would make
    the next test write over the previous one's output.
    """
    def call(src, selector, **opts):
        mod.OPTS.clear()
        mod.OPTS.update({k.replace("_", "-"): str(v) for k, v in opts.items()})
        try:
            return mod.mode_protonate({}, [str(src), selector])
        finally:
            mod.OPTS.clear()
    return call


def pdb_line(serial, name, resname, chain, resseq, xyz, record="ATOM  "):
    """A PDB line with the element columns (77-78) filled, as real files have.

    Without them PDBFile guesses the element from the atom name, so a test built
    on lines that lack them cannot see an element column written wrong - which
    is how a hydrogen once went out declared an oxygen."""
    body = ("%s%5d %-4s%4s %s%4s    %8.3f%8.3f%8.3f  1.00  0.00"
            % (record, serial, name, resname, chain, resseq,
               xyz[0], xyz[1], xyz[2]))
    element = "ZN" if name == "ZN" else name[0]
    return body.ljust(76) + "%2s" % element + "\n"


ASP_SIDECHAIN = [("CB", (0.700, -0.200, 1.000)), ("CG", CG),
                 ("OD1", OD1), ("OD2", OD2)]


@pytest.fixture
def aspartate(tmp_path):
    """An aspartate with a real carboxylate, plus a box, a metal and a water.

    The spectators are here so the pass-through has something to preserve: the
    whole point is that nothing but the one residue moves.
    """
    lines = ["CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1"
             "           1\n"]
    serial = 1
    for name, xyz in [("N", (-1.300, -1.100, 0.000)),
                      ("CA", (0.100, -1.200, 0.000))] + ASP_SIDECHAIN + [
                      ("C", (0.600, -2.600, 0.000)),
                      ("O", (1.800, -2.800, 0.000))]:
        lines.append(pdb_line(serial, name, "ASP", "A", "92", xyz))
        serial += 1
    lines.append(pdb_line(serial, "CA", "GLY", "A", "93", (0.0, -5.0, 0.0)))
    lines.append(pdb_line(serial + 1, "ZN", "ZN", "B", "1", (2.0, 3.0, 0.5),
                          record="HETATM"))
    lines.append(pdb_line(serial + 2, "O", "HOH", "W", "500", (9.0, 9.0, 9.0),
                          record="HETATM"))
    lines.append("END\n")
    path = tmp_path / "start.pdb"
    path.write_text("".join(lines))
    return path


def atoms_of(path):
    """{(chain, resseq, name): (resname, (x, y, z))} over the ATOM records."""
    out = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            out[(line[21], line[22:26].strip(), line[12:16].strip())] = (
                line[17:20].strip(),
                (float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return out


def written(path):
    """The default output path beside `path`."""
    return Path(str(path).rsplit(".", 1)[0] + "_A92_ASH.pdb")


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    return tuple(x / norm(a) for x in a)


# ---- the geometry -----------------------------------------------------------

# Everything read back out of a PDB has been rounded to three decimal places,
# which moves a 0.98 A bond by up to about 0.001 A and the angle by a few
# hundredths of a degree. These tolerances are that rounding and nothing else:
# `exact` below asserts the construction itself to machine precision, so slack
# here cannot hide a wrong placement.
PDB_LENGTH = 2e-3       # A
PDB_ANGLE = 0.2         # degrees


@pytest.fixture
def exact(mod):
    """HD2's position straight from the construction, unrounded."""
    return mod.place_carboxyl_proton(CG, OD1, OD2)


@pytest.fixture
def protonated(run, aspartate):
    """HD2's position as the written PDB actually holds it."""
    run(aspartate, "A:92", to="ASH")
    return atoms_of(written(aspartate))[("A", "92", "HD2")][1]


def test_the_proton_is_an_O_H_bond_length_from_OD2(exact, protonated):
    assert norm(sub(exact, OD2)) == pytest.approx(0.98, abs=1e-12)
    assert norm(sub(protonated, OD2)) == pytest.approx(0.98, abs=PDB_LENGTH)


def test_the_proton_is_in_the_carboxylate_plane(exact, protonated):
    """Out-of-plane distance, measured against the CG/OD1/OD2 normal.

    The carboxylate is sp2 and the proton is conjugated into it; a proton out of
    the plane is a torsion the force field immediately pulls back, which is a
    coordinate the comparison did not mean to perturb."""
    normal = unit(cross(sub(CG, OD2), sub(OD1, OD2)))
    assert abs(dot(sub(exact, OD2), normal)) == pytest.approx(0.0, abs=1e-12)
    assert abs(dot(sub(protonated, OD2), normal)) == pytest.approx(
        0.0, abs=PDB_LENGTH)


def test_the_C_O_H_angle(exact, protonated):
    measure = lambda p: math.degrees(math.acos(
        dot(unit(sub(CG, OD2)), unit(sub(p, OD2)))))
    assert measure(exact) == pytest.approx(113.0, abs=1e-9)
    assert measure(protonated) == pytest.approx(113.0, abs=PDB_ANGLE)


def test_the_proton_is_anti_to_OD1_not_syn(protonated):
    """Anti, asserted as the geometry rather than as a distance to remember.

    Reflecting the built proton in the CG-OD2 axis gives the syn conformer with
    exactly the same bond length and the same C-O-H angle - the only other
    position those two constraints allow in the plane - so the two are directly
    comparable. Anti is the one further from the carbonyl oxygen, and syn is the
    one that would arrive as a 1,3 contact for a minimisation to undo.
    """
    axis = unit(sub(CG, OD2))
    offset = sub(protonated, OD2)
    along = tuple(dot(offset, axis) * a for a in axis)
    perpendicular = sub(offset, along)
    syn = tuple(OD2[i] + along[i] - perpendicular[i] for i in range(3))

    # Same bond length and same angle, so the comparison is between conformers.
    assert norm(sub(syn, OD2)) == pytest.approx(norm(offset), abs=1e-9)
    assert norm(sub(protonated, OD1)) > norm(sub(syn, OD1))

    # And the condition itself: the proton is on the far side of the CG-OD2 axis
    # from OD1, which is what "anti" means.
    to_od1 = sub(OD1, OD2)
    in_plane = sub(to_od1, tuple(dot(to_od1, axis) * a for a in axis))
    assert dot(perpendicular, in_plane) < 0.0


def test_the_placement_follows_the_residue_rather_than_being_written_down(mod):
    """The same carboxylate rotated and translated gives the rotated proton.

    This is what stops a coordinate being hard-coded: the construction is in the
    residue's own frame, so Asp92 wherever it sits in whatever frame of the box
    gets the same proton relative to its own carboxylate.
    """
    here = mod.place_carboxyl_proton(CG, OD1, OD2)
    # A 90 degree turn about z, then a shift.
    turn = lambda p: (-p[1] + 7.0, p[0] - 3.0, p[2] + 11.0)
    there = mod.place_carboxyl_proton(turn(CG), turn(OD1), turn(OD2))
    for got, want in zip(there, turn(here)):
        assert got == pytest.approx(want, abs=1e-9)


def test_a_collinear_carboxylate_raises(mod):
    """No plane, so no in-plane direction and nothing to build along."""
    with pytest.raises(ValueError, match="collinear"):
        mod.place_carboxyl_proton((0.0, 0.0, 0.0), (2.5, 0.0, 0.0),
                                  (1.25, 0.0, 0.0))


# ---- the file ---------------------------------------------------------------

def test_the_residue_is_renamed_and_only_that_residue(run, aspartate):
    run(aspartate, "A:92", to="ASH")
    after = atoms_of(written(aspartate))
    renamed = {key: value[0] for key, value in after.items() if key[1] == "92"}
    assert set(renamed.values()) == {"ASH"}, renamed
    assert after[("A", "93", "CA")][0] == "GLY"
    assert after[("B", "1", "ZN")][0] == "ZN"
    assert after[("W", "500", "O")][0] == "HOH"


def test_every_other_atom_is_copied_through_byte_for_byte(run, aspartate):
    """Including the serial numbers: the proton takes the next free serial
    rather than being inserted into the numbering, so no other line changes and
    no existing CONECT record would have to be rewritten."""
    before = Path(aspartate).read_text().splitlines()
    run(aspartate, "A:92", to="ASH")
    after = written(aspartate).read_text().splitlines()

    untouched = lambda lines: [l for l in lines
                               if not (l.startswith(("ATOM", "HETATM"))
                                       and l[22:26].strip() == "92")
                               and not l.startswith("CONECT")]
    assert untouched(after) == untouched(before)


def test_nothing_outside_the_residue_gains_or_loses_an_atom(run, aspartate):
    before, _ = atoms_of(aspartate), run(aspartate, "A:92", to="ASH")
    after = atoms_of(written(aspartate))
    elsewhere = lambda d: {k for k in d if k[1] != "92"}
    assert elsewhere(after) == elsewhere(before)
    assert len(after) == len(before) + 1
    assert set(after) - set(before) == {("A", "92", "HD2")}


def test_the_box_record_survives(run, aspartate):
    run(aspartate, "A:92", to="ASH")
    assert written(aspartate).read_text().startswith("CRYST1   40.000")


def test_the_proton_is_declared_a_hydrogen_in_the_element_columns(run, aspartate):
    """PDBFile believes columns 77-78 over the atom name. The line is built from
    OD2's, so without this the proton would arrive declared an oxygen and the
    residue would be an aspartate with three sidechain oxygens."""
    run(aspartate, "A:92", to="ASH")
    line = next(l for l in written(aspartate).read_text().splitlines()
                if l[12:16].strip() == "HD2")
    assert line[76:78] == " H"


def test_the_proton_sits_inside_its_residue_and_carries_a_fresh_serial(run,
                                                                      aspartate):
    """PDBFile groups atoms into residues by file order, so a proton appended
    at the end of the file would become a residue of its own."""
    run(aspartate, "A:92", to="ASH")
    lines = [l for l in written(aspartate).read_text().splitlines()
             if l.startswith(("ATOM", "HETATM"))]
    resids = [l[22:26].strip() for l in lines]
    first, last = resids.index("92"), len(resids) - 1 - resids[::-1].index("92")
    assert resids[first:last + 1] == ["92"] * (last + 1 - first)

    serials = [int(l[6:11]) for l in lines]
    assert len(set(serials)) == len(serials), "a duplicate serial breaks CONECT"
    proton = next(l for l in lines if l[12:16].strip() == "HD2")
    assert int(proton[6:11]) == max(int(l[6:11]) for l in lines
                                    if l[12:16].strip() != "HD2") + 1


def test_a_CONECT_record_states_the_new_bond_before_END(run, aspartate):
    """A CONECT after END is dropped without complaint, and an unbonded proton
    matches no template."""
    run(aspartate, "A:92", to="ASH")
    lines = written(aspartate).read_text().splitlines()
    serial = {l[12:16].strip(): int(l[6:11]) for l in lines
              if l.startswith(("ATOM", "HETATM")) and l[22:26].strip() == "92"}
    record = "CONECT%5d%5d" % (serial["OD2"], serial["HD2"])
    assert record in lines
    assert lines.index(record) < lines.index("END")


def test_the_input_is_never_the_output(run, aspartate, tmp_path):
    before = Path(aspartate).read_text()
    with pytest.raises(SystemExit, match="does not edit one in place"):
        run(aspartate, "A:92", to="ASH", out=str(aspartate))
    assert Path(aspartate).read_text() == before


def test_out_is_honoured(run, aspartate, tmp_path):
    out = tmp_path / "elsewhere" / "ash.pdb"
    out.parent.mkdir()
    run(aspartate, "A:92", to="ASH", out=str(out))
    assert ("A", "92", "HD2") in atoms_of(out)
    assert not written(aspartate).exists()


# ---- the refusals -----------------------------------------------------------

def test_a_residue_that_is_not_an_aspartate_raises(run, aspartate, tmp_path):
    src = tmp_path / "asn.pdb"
    src.write_text(Path(aspartate).read_text().replace("ASP A  92",
                                                       "ASN A  92"))
    with pytest.raises(SystemExit, match="ASN.*ASH is built from ASP"):
        run(src, "A:92", to="ASH")
    assert not written(src).exists()


def test_a_residue_that_already_has_HD2_raises(run, aspartate, tmp_path):
    """Running it twice must not build a second proton on the same oxygen."""
    src = tmp_path / "already.pdb"
    lines = Path(aspartate).read_text().splitlines(keepends=True)
    at = next(i for i, l in enumerate(lines) if l[12:16].strip() == "OD2")
    lines.insert(at + 1, pdb_line(99, "HD2", "ASP", "A", "92",
                                  (0.186, -0.666, 0.800)))
    src.write_text("".join(lines))
    with pytest.raises(SystemExit, match="already carries HD2"):
        run(src, "A:92", to="ASH")
    assert not written(src).exists()


def test_the_output_cannot_be_protonated_again(run, aspartate):
    """The residue is ASH by then, which is not what ASH is built from."""
    run(aspartate, "A:92", to="ASH")
    with pytest.raises(SystemExit, match="ASH is built from ASP"):
        run(written(aspartate), "A:92", to="ASH")


def test_a_missing_residue_raises(run, aspartate):
    with pytest.raises(SystemExit, match="no residue A:999"):
        run(aspartate, "A:999", to="ASH")


def test_a_carboxylate_missing_an_oxygen_raises(run, aspartate, tmp_path):
    src = tmp_path / "broken.pdb"
    src.write_text("".join(l for l in Path(aspartate).read_text()
                           .splitlines(keepends=True)
                           if l[12:16].strip() != "OD1"))
    with pytest.raises(SystemExit, match="missing OD1"):
        run(src, "A:92", to="ASH")


def test_an_unknown_variant_raises(run, aspartate):
    with pytest.raises(SystemExit, match="--to="):
        run(aspartate, "A:92", to="ASX")
    with pytest.raises(SystemExit, match="--to="):
        run(aspartate, "A:92")


def test_LYN_is_refused_with_the_reason(run, aspartate):
    """Not half-done: ASH and GLH work by renaming because pdbNames.xml maps
    them back to ASP/GLU, whose residues.xml entries already carry HD2 and HE2.
    Nothing maps LYN, and residues.xml has no LYN entry, so a residue renamed
    LYN comes out of createStandardBonds with no bonds at all. Neutralising a
    lysine is `make_water_nuc.py --lys`, which leaves the name LYS."""
    with pytest.raises(SystemExit, match="not implemented"):
        run(aspartate, "A:92", to="LYN")


def test_a_bad_selector_raises(run, aspartate):
    with pytest.raises(SystemExit, match="CHAIN:RESID"):
        run(aspartate, "92", to="ASH")


# ---- against the force field ------------------------------------------------
#
# The point of the whole mode: a residue that ff14SB can template. Asserted
# against the force field rather than against a list written out here, because a
# list in a test is just the same assumption twice.

CAPS = [("HH31", "ACE", "1", (-4.0, 1.0, 0.0)), ("CH3", "ACE", "1", (-3.5, 0.0, 0.0)),
        ("HH32", "ACE", "1", (-4.0, -1.0, 0.0)), ("HH33", "ACE", "1", (-3.5, 0.0, 1.0)),
        ("C", "ACE", "1", (-2.0, 0.0, 0.0)), ("O", "ACE", "1", (-1.5, 1.1, 0.0))]
TAIL = [("N", "NME", "93", (-0.3, -3.6, 0.0)), ("H", "NME", "93", (-1.3, -3.4, 0.0)),
        ("CH3", "NME", "93", (0.1, -5.0, 0.0)), ("HH31", "NME", "93", (1.2, -5.1, 0.0)),
        ("HH32", "NME", "93", (-0.3, -5.5, 0.9)),
        ("HH33", "NME", "93", (-0.3, -5.5, -0.9))]
BODY = [("N", (-1.300, -1.100, 0.000)), ("H", (-1.800, -2.000, 0.000)),
        ("CA", (0.100, -1.200, 0.000)), ("HA", (0.400, -2.200, 0.300)),
        ("CB", (0.700, -0.200, 1.000)), ("HB2", (0.400, -0.500, 2.000)),
        ("HB3", (1.800, -0.200, 1.000)), ("CG", CG), ("OD1", OD1), ("OD2", OD2),
        ("C", (0.600, -2.600, 0.000)), ("O", (1.800, -2.800, 0.000))]


@pytest.fixture
def capped_aspartate(tmp_path):
    """ACE-ASP-NME: ff14SB's ASH is a mid-chain template, with no N- or
    C-terminal variant, so the residue needs a neighbour at each end to be
    matched at all. The backbone coordinates are crude on purpose -
    `createStandardBonds` bonds by atom name with no distance test, and nothing
    here is minimised - while the carboxylate is the real geometry."""
    rows = CAPS + [(name, "ASP", "92", xyz) for name, xyz in BODY] + TAIL
    lines = [pdb_line(i, name, resname, "A", resseq, xyz)
             for i, (name, resname, resseq, xyz) in enumerate(rows, 1)]
    path = tmp_path / "capped.pdb"
    path.write_text("".join(lines) + "END\n")
    return path


def test_the_atom_names_are_exactly_the_ff14SB_ASH_template(run,
                                                            capped_aspartate):
    pytest.importorskip("openmm", reason="the template lives in ff14SB's XML")
    from openmm.app import ForceField, PDBFile

    run(capped_aspartate, "A:92", to="ASH")
    out = Path(str(capped_aspartate).rsplit(".", 1)[0] + "_A92_ASH.pdb")
    field = ForceField("amber14/protein.ff14SB.xml")
    pdb = PDBFile(str(out))

    residue = next(r for r in pdb.topology.residues() if r.id.strip() == "92")
    # PDBFile renames it: pdbNames.xml lists ASH as an alternative name for ASP,
    # so the topology says ASP and the ASH *template* is matched afterwards by
    # bond graph. That is also why `amino_acids` never sees the name ASH.
    assert residue.name == "ASP"
    names = [atom.name for atom in residue.atoms()]
    template = [atom.name for atom in field._templates["ASH"].atoms]
    assert set(names) == set(template)
    assert len(names) == len(template)

    # And the proton is bonded, which is what the template is matched on.
    proton = next(a for a in residue.atoms() if a.name == "HD2")
    partners = [b[0].name if b[1] is proton else b[1].name
                for b in pdb.topology.bonds() if proton in b]
    assert partners == ["OD2"]


def test_the_structure_builds_a_system_and_the_residue_is_neutral(
        run, capped_aspartate):
    """The failure this mode exists to prevent is `createSystem` raising "No
    template found for residue" hours into a job, so the test is that call."""
    pytest.importorskip("openmm", reason="createSystem is OpenMM's")
    import openmm
    from openmm.app import ForceField, PDBFile

    field = ForceField("amber14/protein.ff14SB.xml")
    run(capped_aspartate, "A:92", to="ASH")
    out = Path(str(capped_aspartate).rsplit(".", 1)[0] + "_A92_ASH.pdb")

    charge_of = {}
    for path, label in ((capped_aspartate, "ASP"), (out, "ASH")):
        pdb = PDBFile(str(path))
        system = field.createSystem(pdb.topology)
        nonbonded = next(f for f in system.getForces()
                         if isinstance(f, openmm.NonbondedForce))
        residue = next(r for r in pdb.topology.residues()
                       if r.id.strip() == "92")
        charge_of[label] = sum(
            nonbonded.getParticleParameters(a.index)[0]._value
            for a in residue.atoms())

    assert charge_of["ASP"] == pytest.approx(-1.0, abs=1e-6)
    assert charge_of["ASH"] == pytest.approx(0.0, abs=1e-6)


# ---- GLH, the same construction one atom further out ------------------------

GLU_CD = (0.800, 2.000, -0.300)
OE1 = (GLU_CD[0] + 1.25 * math.cos(HALF_OCO),
       GLU_CD[1] + 1.25 * math.sin(HALF_OCO), GLU_CD[2])
OE2 = (GLU_CD[0] + 1.25 * math.cos(-HALF_OCO),
       GLU_CD[1] + 1.25 * math.sin(-HALF_OCO), GLU_CD[2])

GLU_BODY = [("N", (-1.300, -1.100, 0.000)), ("H", (-1.800, -2.000, 0.000)),
            ("CA", (0.100, -1.200, 0.000)), ("HA", (0.400, -2.200, 0.300)),
            ("CB", (0.700, -0.200, 1.000)), ("HB2", (0.400, -0.500, 2.000)),
            ("HB3", (1.800, -0.200, 1.000)), ("CG", (0.200, 1.200, 0.900)),
            ("HG2", (-0.900, 1.200, 0.900)), ("HG3", (0.500, 1.700, 1.800)),
            ("CD", GLU_CD), ("OE1", OE1), ("OE2", OE2),
            ("C", (0.600, -2.600, 0.000)), ("O", (1.800, -2.800, 0.000))]


@pytest.fixture
def capped_glutamate(tmp_path):
    rows = ([(n, r, s, x) for n, r, s, x in CAPS]
            + [(name, "GLU", "150", xyz) for name, xyz in GLU_BODY]
            + [(n, r, "151", x) for n, r, _, x in TAIL])
    lines = [pdb_line(i, name, resname, "A", resseq, xyz)
             for i, (name, resname, resseq, xyz) in enumerate(rows, 1)]
    path = tmp_path / "glu.pdb"
    path.write_text("".join(lines) + "END\n")
    return path


def test_GLH_goes_the_same_way_on_OE2(run, capped_glutamate):
    """The entry is a table of atom names, so GLH is the same code: HE2 on OE2,
    built from CD/OE1/OE2, and the same template check applies."""
    pytest.importorskip("openmm", reason="the template lives in ff14SB's XML")
    from openmm.app import ForceField, PDBFile

    run(capped_glutamate, "A:150", to="GLH")
    out = Path(str(capped_glutamate).rsplit(".", 1)[0] + "_A150_GLH.pdb")
    proton = atoms_of(out)[("A", "150", "HE2")][1]
    assert norm(sub(proton, OE2)) == pytest.approx(0.98, abs=PDB_LENGTH)
    assert math.degrees(math.acos(dot(unit(sub(GLU_CD, OE2)),
                                      unit(sub(proton, OE2))))) \
        == pytest.approx(113.0, abs=PDB_ANGLE)
    normal = unit(cross(sub(GLU_CD, OE2), sub(OE1, OE2)))
    assert abs(dot(sub(proton, OE2), normal)) == pytest.approx(
        0.0, abs=PDB_LENGTH)
    assert norm(sub(proton, OE1)) > norm(sub(OE2, OE1))

    field = ForceField("amber14/protein.ff14SB.xml")
    pdb = PDBFile(str(out))
    residue = next(r for r in pdb.topology.residues() if r.id.strip() == "150")
    assert {a.name for a in residue.atoms()} == {
        a.name for a in field._templates["GLH"].atoms}
    field.createSystem(pdb.topology)
