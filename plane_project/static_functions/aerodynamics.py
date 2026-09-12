import aerosandbox.numpy as np
import aerosandbox as asb
from aerosandbox.tools import units


def create_wing(v):
    def const_quarter_chord(rc, tc):
        return rc / 4 - tc / 4

    tip_z = v.b / 2 * np.sin(v.dihedral * np.pi / 180)

    v.wing = asb.Wing(
        name="wing",
        symmetric=True,
        xsecs=[
            asb.WingXSec(
                xyz_le=[0, 0, 0],  # UPDATED per Taras/AVL documentation
                chord=v.c_r,
                twist=0,
                airfoil=v.airfoil,
                control_surfaces=[
                    asb.ControlSurface(
                        name="Aileron",
                        hinge_point=.7,
                        deflection=0
                    )
                ]
            ),
            asb.WingXSec(
                xyz_le=[
                    const_quarter_chord(v.c_r, v.c_t),
                    v.b / 2,
                    tip_z
                ],
                chord=v.c_t,
                twist=v.washout,
                airfoil=v.airfoil
            )
        ]
    )


def create_hor_stab(v):
    v.hor_stab = asb.Wing(
        name="horizontal_tail",
        symmetric=True,
        xsecs=[
            asb.WingXSec(
                xyz_le=[0, 0, 0],
                chord=v.c_h,
                twist=v.decalage,
                airfoil=v.airfoil_tail,
                control_surfaces=[
                    asb.ControlSurface(
                        name="elevator",
                        hinge_point=0.70,
                        deflection=0
                    )
                ]
            ),
            asb.WingXSec(
                xyz_le=[0, v.b_h / 2, 0],
                chord=v.c_h,
                twist=v.decalage,
                airfoil=v.airfoil_tail
            )
        ]
    ).translate([v.x_tail, 0, 0])


def create_vert_stab(v):
    v.vert_stab = asb.Wing(
        name="vertical tail",
        xsecs=[
            asb.WingXSec(
                xyz_le=[0, 0, 0],
                chord=v.c_v,
                twist=0,
                airfoil=v.airfoil_tail,
                control_surfaces=[
                    asb.ControlSurface(
                        name="rudder",
                        hinge_point=0.75,
                        deflection=0
                    )
                ]
            ),
            asb.WingXSec(
                xyz_le=[0, 0, v.b_v],
                chord=v.c_v,
                twist=0,
                airfoil=v.airfoil_tail
            )
        ]
    ).translate([v.x_tail, 0, 0])

def run_aero_constraints(v):
    constraints_tail_dim(v)
    constraints_wing_dim(v)

def constraints_wing_dim(v):
    v.opti.subject_to(v.b < 5.8 * units.foot)


def constraints_tail_dim(v):
    v.opti.subject_to(v.c_v == v.c_h)
    v.opti.subject_to(v.x_tail > v.c_r)


def op_point(v):
    v.velocity_M2 = v.opti.variable(
        init_guess=20,
        lower_bound=0,
        upper_bound=85 * units.mph
    )

    v.velocity_M3 = v.opti.variable(
        init_guess=20,
        lower_bound=0,
        upper_bound=85 * units.mph
    )

    v.q_M2 = 1 / 2 * v.rho * v.velocity_M2**2
    v.q_M3 = 1 / 2 * v.rho * v.velocity_M3**2

    v.alpha_M2 = v.opti.variable(
        init_guess=0,
        upper_bound=13,
        lower_bound=-10
    )

    v.alpha_M3 = v.opti.variable(
        init_guess=0,
        upper_bound=13,
        lower_bound=-10
    )

    v.op_pt_M2 = asb.OperatingPoint(
        velocity=v.velocity_M2,
        alpha=v.alpha_M2,
        beta=0
    )

    v.op_pt_M3 = asb.OperatingPoint(
        velocity=v.velocity_M3,
        alpha=v.alpha_M3,
        beta=0
    )


def eval_aero(v):
    # returns an AeroBuildup
    v.plane = asb.Airplane(
        name="Prototype 0",
        xyz_ref=[v.r_c / 4, 0, 0],
        wings=[v.wing, v.hor_stab, v.vert_stab]
    )

    v.ab2 = asb.AeroBuildup(
        airplane=v.plane,
        op_point=v.op_pt_M2
    )

    v.ab3 = asb.AeroBuildup(
        airplane=v.plane,
        op_point=v.op_pt_M3
    )

    v.aero_M2 = v.ab2.run_with_stability_derivatives(
        alpha=True,
        beta=True,
        p=False,
        q=False,
        r=False
    )

    v.aero_M3 = v.ab3.run_with_stability_derivatives(
        alpha=True,
        beta=True,
        p=False,
        q=False,
        r=False
    )


def trim_conditions(opti, params, aero_buildup):
    # finds the trim of the elevator for a given case
    # and maybe aileron deflection depending on turns
    # start a new subopti just to find the trim?
    # return bad if not applicable?
    pass


def passive_stability(opti, params, aero_buildup):
    # checks the passive stability

    # takes the AeroBuildup and runs the derivatives on them,
    # including corrections for the fuselage
    pass


def control_authority(opti, params):
    # check turning rate of ailerons and stuff like that?
    # check how fast it can turn the vehicle in cases of wind?
    # check wing loading and wind effects
    # --> what is an acceptable amount to get blown around by the wind?
    pass


def takeoff_conditions(opti, params):
    # check takeoff length,
    # make sure takeoff is achievable i think?
    pass