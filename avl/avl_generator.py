"""
Generate an AVL geometry file.

Change the parameters in the GEOMETRY section, then run:

    python make_avl.py

It will create:
    prototype_0.avl
"""

from pathlib import Path
import math


# =============================================================================
# GEOMETRY
# =============================================================================

# -------------------------------------------------------------------------
# Wing
# -------------------------------------------------------------------------

WING_SPAN = 1.76784          # full span [m]
WING_ROOT_CHORD = 0.48      # root chord [m]
WING_TAPER = 1.0             # tip chord / root chord
WING_SWEEP = 0.0             # quarter-chord sweep [deg]
WING_DIHEDRAL = 0.0          # dihedral [deg]
WING_INCIDENCE = 0.0         # incidence [deg]

WING_AIRFOIL = "sd7032.dat"

# AVL discretization
WING_NCHORD = 12
WING_NSPAN = 12
WING_CSPACE = 1
WING_SSPACE = 1

# Flap
FLAP_HINGE = 0.60

# Aileron
AILERON_HINGE = 0.70


# -------------------------------------------------------------------------
# Horizontal tail
# -------------------------------------------------------------------------

HT_SPAN = 0.50286508         # full span [m]
HT_ROOT_CHORD = 0.16020765
HT_TAPER = 1.0
HT_SWEEP = 0.0
HT_DIHEDRAL = 0.0
HT_INCIDENCE = -5.0

HT_XLE = 1.0
HT_ZLE = 0.25143254

HT_AIRFOIL = "naca0012.dat"

HT_NCHORD = 12
HT_NSPAN = 12
HT_CSPACE = 1
HT_SSPACE = 1

ELEVATOR_HINGE = 0.70


# -------------------------------------------------------------------------
# Vertical tail
# -------------------------------------------------------------------------

VT_HEIGHT = 0.25143254       # vertical span [m]
VT_ROOT_CHORD = 0.16020765
VT_TAPER = 1.0
VT_SWEEP = 0.0
VT_INCIDENCE = 0.0

VT_XLE = 1.0
VT_AIRFOIL = "naca0012.dat"

VT_NCHORD = 12
VT_NSPAN = 12
VT_CSPACE = 1
VT_SSPACE = 1

RUDDER_HINGE = 0.75


# -------------------------------------------------------------------------
# Reference geometry
# -------------------------------------------------------------------------

SREF = 0.8548320475195411
CREF = 0.4835460491444603
BREF = 1.767839999999997

XREF = 0.0
YREF = 0.0
ZREF = 0.0

MACH = 0.0


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def format_number(x):
    """Format a number compactly for AVL."""
    return f"{x:.10g}"


def wing_tip_chord(root_chord, taper):
    return root_chord * taper


def quarter_chord_x(
    y,
    root_chord,
    taper,
    span,
    sweep_deg,
):
    """
    Calculate x-location of the quarter-chord line.

    y is the full-span coordinate from the aircraft centerline.
    """

    half_span = span / 2

    sweep_rad = math.radians(sweep_deg)

    # Quarter-chord x displacement
    x_qc = abs(y) * math.tan(sweep_rad)

    return x_qc


def leading_edge_x(
    y,
    root_chord,
    taper,
    span,
    sweep_deg,
):
    """
    Calculate leading-edge x coordinate from quarter-chord sweep.
    """

    c = root_chord + (
        root_chord * taper - root_chord
    ) * abs(y) / (span / 2)

    x_qc = quarter_chord_x(
        y,
        root_chord,
        taper,
        span,
        sweep_deg,
    )

    return x_qc - 0.25 * root_chord + 0.25 * c


def wing_section(
    y,
    root_chord,
    taper,
    span,
    sweep_deg,
    dihedral_deg,
    incidence_deg,
):
    """
    Return AVL SECTION geometry at span location y.
    """

    half_span = span / 2

    # Chord varies linearly with span
    eta = abs(y) / half_span

    chord = root_chord * (
        1 + (taper - 1) * eta
    )

    # Leading edge x
    x = leading_edge_x(
        y,
        root_chord,
        taper,
        span,
        sweep_deg,
    )

    # Dihedral
    z = abs(y) * math.tan(
        math.radians(dihedral_deg)
    )

    return x, y, z, chord, incidence_deg


def add_section(
    lines,
    section,
    airfoil,
    claf,
    control_lines=None,
):
    """
    Add one AVL SECTION block.
    """

    x, y, z, chord, incidence = section

    lines.append("#--------------------------------------------------")
    lines.append("SECTION")
    lines.append(
        "#Xle    Yle    Zle     Chord   Ainc"
    )
    lines.append(
        f"{format_number(x)} "
        f"{format_number(y)} "
        f"{format_number(z)} "
        f"{format_number(chord)} "
        f"{format_number(incidence)}"
    )
    lines.append("")
    lines.append("AFIL")
    lines.append(airfoil)
    lines.append("")
    lines.append("CLAF")
    lines.append(
        f"{format_number(claf)}"
        "  # Computed using rule from avl_doc.txt"
    )

    if control_lines:
        for control in control_lines:
            lines.append("")
            lines.append("CONTROL")
            lines.append(
                "#name, gain, Xhinge, XYZhvec, SgnDup"
            )
            lines.append(control)


def generate_wing(lines):

    root_section = wing_section(
        y=0,
        root_chord=WING_ROOT_CHORD,
        taper=WING_TAPER,
        span=WING_SPAN,
        sweep_deg=WING_SWEEP,
        dihedral_deg=WING_DIHEDRAL,
        incidence_deg=WING_INCIDENCE,
    )

    mid_section = wing_section(
        y=WING_SPAN/4,
        root_chord=WING_ROOT_CHORD,
        taper=WING_TAPER,
        span=WING_SPAN,
        sweep_deg=WING_SWEEP,
        dihedral_deg=WING_DIHEDRAL,
        incidence_deg=WING_INCIDENCE,
    )


    tip_section = wing_section(
        y=WING_SPAN / 2,
        root_chord=WING_ROOT_CHORD,
        taper=WING_TAPER,
        span=WING_SPAN,
        sweep_deg=WING_SWEEP,
        dihedral_deg=WING_DIHEDRAL,
        incidence_deg=WING_INCIDENCE,
    )

    lines.append("#===============================================================================")
    lines.append("SURFACE")
    lines.append("wing")
    lines.append("#Nchordwise  Cspace  [Nspanwise   Sspace]")
    lines.append(
        f"{WING_NCHORD} {WING_CSPACE} "
        f"{WING_NSPAN} {WING_SSPACE}"
    )
    lines.append("")
    lines.append("YDUPLICATE")
    lines.append("0")
    lines.append("")

    # Root
    add_section(
        lines,
        root_section,
        WING_AIRFOIL,
        claf=1.076655270298305,
        control_lines=[
            f"flap 1 {FLAP_HINGE} 0 0 0 1"
        ],
    )
    add_section(
        lines,
        mid_section,
        WING_AIRFOIL,
        claf=1.076655270298305,
        control_lines=[
            f"flap 1 {FLAP_HINGE} 0 0 0 1\n"
            f"aileron 1 {AILERON_HINGE} 0 0 0 -1",
        ],
    )

    # Tip
    add_section(
        lines,
        tip_section,
        WING_AIRFOIL,
        claf=1.076655270298305,
        control_lines=[
            f"flap 1 {FLAP_HINGE} 0 0 0 1\n",
            f"aileron 1 {AILERON_HINGE} 0 0 0 -1",
        ],
    )


def generate_horizontal_tail(lines):

    root = [
        HT_XLE,
        0,
        HT_ZLE,
        HT_ROOT_CHORD,
        HT_INCIDENCE,
    ]

    tip_chord = HT_ROOT_CHORD * HT_TAPER

    tip_x = (
        HT_XLE
        + (HT_SPAN / 2)
        * math.tan(math.radians(HT_SWEEP))
    )

    tip_z = (
        HT_ZLE
        + (HT_SPAN / 2)
        * math.tan(math.radians(HT_DIHEDRAL))
    )

    tip = [
        tip_x,
        HT_SPAN / 2,
        tip_z,
        tip_chord,
        HT_INCIDENCE,
    ]

    lines.append("#===============================================================================")
    lines.append("SURFACE")
    lines.append("horizontal_tail")
    lines.append("#Nchordwise  Cspace  [Nspanwise   Sspace]")
    lines.append(
        f"{HT_NCHORD} {HT_CSPACE} "
        f"{HT_NSPAN} {HT_SSPACE}"
    )
    lines.append("")
    lines.append("YDUPLICATE")
    lines.append("0")
    lines.append("")

    add_section(
        lines,
        root,
        HT_AIRFOIL,
        claf=1.0924221254554969,
        control_lines=[
            f"elevator 1 {ELEVATOR_HINGE} 0 0 0 1"
        ],
    )

    add_section(
        lines,
        tip,
        HT_AIRFOIL,
        claf=1.0924221254554969,
        control_lines=[
            f"elevator 1 {ELEVATOR_HINGE} 0 0 0 1"
        ],
    )


def generate_vertical_tail(lines):

    tip_chord = VT_ROOT_CHORD * VT_TAPER

    tip_x = (
        VT_XLE
        + VT_HEIGHT
        * math.tan(math.radians(VT_SWEEP))
    )

    root = [
        VT_XLE,
        0,
        0,
        VT_ROOT_CHORD,
        VT_INCIDENCE,
    ]

    tip = [
        tip_x,
        0,
        VT_HEIGHT,
        tip_chord,
        VT_INCIDENCE,
    ]

    lines.append("#===============================================================================")
    lines.append("SURFACE")
    lines.append("vertical tail")
    lines.append("#Nchordwise  Cspace  [Nspanwise   Sspace]")
    lines.append(
        f"{VT_NCHORD} {VT_CSPACE} "
        f"{VT_NSPAN} {VT_SSPACE}"
    )
    lines.append("")
    lines.append("")

    add_section(
        lines,
        root,
        VT_AIRFOIL,
        claf=1.0924221254554969,
        control_lines=[
            f"rudder 1 {RUDDER_HINGE} 0 0 0 1"
        ],
    )

    add_section(
        lines,
        tip,
        VT_AIRFOIL,
        claf=1.0924221254554969,
        control_lines=[
            f"rudder 1 {RUDDER_HINGE} 0 0 0 1"
        ],
    )


# =============================================================================
# WRITE AVL FILE
# =============================================================================

def generate_avl():

    lines = []

    lines.append("Prototype 0")
    lines.append("#Mach")
    lines.append(
        f"{format_number(MACH)}"
        "        ! AeroSandbox note: This is overwritten later "
        "to match the current OperatingPoint Mach during the AVL run."
    )

    lines.append("#IYsym   IZsym   Zsym")
    lines.append("0       0   0")

    lines.append("#Sref    Cref    Bref")
    lines.append(
        f"{format_number(SREF)} "
        f"{format_number(CREF)} "
        f"{format_number(BREF)}"
    )

    lines.append("#Xref    Yref    Zref")
    lines.append(
        f"{format_number(XREF)} "
        f"{format_number(YREF)} "
        f"{format_number(ZREF)}"
    )

    lines.append("# CDp")
    lines.append("0")
    lines.append("")

    generate_wing(lines)
    generate_horizontal_tail(lines)
    generate_vertical_tail(lines)

    # ---------------------------------------------------------------------
    # Optional fuselage
    # ---------------------------------------------------------------------

    lines.append("")
    lines.append("#==================================================")
    lines.append("#BODY")
    lines.append("#Untitled")
    lines.append("#24 1")
    lines.append("")
    lines.append("#BFIL")
    lines.append("#M2.avl.fuse1")
    lines.append("")
    lines.append("#TRANSLATE")
    lines.append("#0 0 0")

    return "\n".join(lines) + "\n"


if __name__ == "__main__":

    output_file = Path("prototype_0.avl")

    avl_text = generate_avl()

    output_file.write_text(avl_text)

    print(f"Wrote {output_file}")
