import aerosandbox as asb
import aerosandbox.numpy as np
import json
from aerosandbox.tools import units
import re
from pathlib import Path
import pandas as pd
import numpy as onp   # plain numpy (np above is aerosandbox.numpy, which is for symbolic math)
import casadi as ca
import time

COLUMNS = ["alpha", "CL", "CD", "CDp", "CM", "Top_Xtr", "Bot_Xtr", "Top_Itr", "Bot_Itr"]
NAME_RE = re.compile(r"^(?P<airfoil>.+)_Re(?P<Re>\d+(?:\.\d+)?)$")  # greedy, so underscores in the airfoil name are OK

def parse_polar_rows(path):
    """Read the numeric rows (after the dashed line) of one XFOIL polar file."""
    lines = Path(path).read_text(errors="ignore").splitlines()
    start = next((i for i, l in enumerate(lines) if l.strip().startswith("------")), None)
    rows = []
    if start is not None:
        for l in lines[start + 1:]:
            parts = l.split()
            if len(parts) == len(COLUMNS):
                try:
                    rows.append([float(x) for x in parts])
                except ValueError:
                    pass
    return pd.DataFrame(rows, columns=COLUMNS)

def load_airfoils(folder, pattern="*.txt"):
    """Returns {airfoil_name: DataFrame}, with all Reynolds numbers stacked in one table."""
    frames = {}
    for p in sorted(Path(folder).glob(pattern)):
        m = NAME_RE.match(p.stem)
        if not m:
            print(f"Skipping (name doesn't match pattern): {p.name}")
            continue
        df = parse_polar_rows(p)
        df.insert(0, "Re", float(m.group("Re")))
        frames.setdefault(m.group("airfoil"), []).append(df)

    return {
        name: pd.concat(dfs, ignore_index=True).sort_values(["Re", "alpha"]).reset_index(drop=True)
        for name, dfs in frames.items()
    }

class AirfoilPolar:
    """
    Smooth, differentiable CL(alpha) and CD(alpha) for ONE airfoil at ONE Reynolds number,
    built from the DataFrame that load_airfoils() returns (columns: Re, alpha, CL, CD, ...).

    - If Re falls between two Reynolds numbers in the data, the polars are blended in log(Re).
    - Only the pre-stall branch is kept (alpha up to CL_max); post-stall XFOIL data is unreliable.
    - Evaluated as a cubic B-spline (casadi), so it works inside asb.Opti and IPOPT gets real
      derivatives. (Plain linear interpolation has kinks and made IPOPT fail in testing.)
    """
    def __init__(self, df, Re, n_grid=60):
        self.Re_requested = float(Re)
        res = onp.sort(df["Re"].unique())
        if len(res) == 0:
            raise ValueError("polar file(s) contain no converged data points")
        self.Re = float(onp.clip(Re, res.min(), res.max()))
        self.Re_clamped = bool(self.Re != Re)  # True if the data didn't cover the requested Re

        lo, hi = res[res <= self.Re].max(), res[res >= self.Re].min()
        w = 0.0 if lo == hi else (onp.log(self.Re) - onp.log(lo)) / (onp.log(hi) - onp.log(lo))

        def clean(r):
            d = df[df.Re == r].sort_values("alpha").drop_duplicates("alpha")
            return d[onp.isfinite(d.CL) & onp.isfinite(d.CD) & (d.CD > 0)]
        d_lo, d_hi = clean(lo), clean(hi)
        if len(d_lo) < 8 or len(d_hi) < 8:
            raise ValueError("too few converged polar points")

        # blend the two bracketing polars on their common alpha range
        amin = max(d_lo.alpha.min(), d_hi.alpha.min())
        amax = min(d_lo.alpha.max(), d_hi.alpha.max())
        a = onp.linspace(amin, amax, 300)
        cl = (1 - w) * onp.interp(a, d_lo.alpha, d_lo.CL) + w * onp.interp(a, d_hi.alpha, d_hi.CL)
        cd = (1 - w) * onp.interp(a, d_lo.alpha, d_lo.CD) + w * onp.interp(a, d_hi.alpha, d_hi.CD)

        # keep the pre-stall branch only, then resample to a uniform grid for the spline
        i_stall = int(onp.argmax(cl))
        self.alpha_min, self.alpha_stall, self.CL_max = float(a[0]), float(a[i_stall]), float(cl[i_stall])
        if i_stall < 10:
            raise ValueError("polar has no usable pre-stall range")
        grid = onp.linspace(self.alpha_min, self.alpha_stall, n_grid)
        self.alpha_grid = grid
        self.CL_grid = onp.interp(grid, a[: i_stall + 1], cl[: i_stall + 1])
        self.CD_grid = onp.interp(grid, a[: i_stall + 1], cd[: i_stall + 1])
        self.CD_min = float(self.CD_grid.min())

        self._cl = ca.interpolant("cl", "bspline", [grid], self.CL_grid)
        self._cd = ca.interpolant("cd", "bspline", [grid], self.CD_grid)

    def cl(self, alpha):  # alpha in degrees; works with opti variables
        return self._cl(alpha)

    def cd(self, alpha):
        return self._cd(alpha)


class variables:
    def __init__(self, airfoil=None, polar_df=None, Re=500_000, stall_margin=0.9):
        #airfoil:      UIUC name (str) or asb.Airfoil -> the real airfoil is used for the weight/structure sizing
        #polar_df:     airfoils[name] from load_airfoils() -> CL and CD come from this polar instead of the placeholder model
        #Re:           Reynolds number to read the polar at. One number for everything, or {"M2": ..., "M3": ...}
        #stall_margin: stay below this fraction of the polar's CL_max in flight
        self.given_airfoil = airfoil is not None
        self.opti = asb.Opti()
        self.order = 2
        self.mission_no=0

        #overall wing design
        self.S = self.opti.variable(init_guess = 0.2, lower_bound = 0.01, upper_bound = 1)
        self.AR = self.opti.variable(init_guess = 3, lower_bound = 2, upper_bound = 15)
        # self.AR = (5.79 * units.foot) ** 2 / self.S
        self.taper = 1
        if airfoil is None:
            airfoil = "sd7032" #default for the no-airfoil path (only used for the aero-plane geometry)
        self.airfoil = airfoil if isinstance(airfoil, asb.Airfoil) else asb.Airfoil(airfoil.lower())
        if self.airfoil.coordinates is None: #asb only warns when a name isn't in the UIUC database
            raise ValueError(f"airfoil '{airfoil}' not found in the aerosandbox UIUC database")

        self.CL = {} #filled per mission by CL_no_af
        self.CD = {} #filled per mission by CD_plane_no_af

        self.dihedral = 0 #degrees
        self.washout = 0 #degrees
        #tail design
        self.AR_v = self.opti.variable(init_guess = 2, lower_bound = 1, upper_bound = 5)
        self.AR_h = self.opti.variable(init_guess = 4, lower_bound = 3, upper_bound = 5)
        self.S_v = self.opti.variable(init_guess = 0.01, lower_bound = 0.001, upper_bound = .1)
        self.S_h = self.opti.variable(init_guess = 0.01, lower_bound = 0.001, upper_bound = .1)

        self.decalage = -5 #degree

        self.airfoil_tail = asb.Airfoil("naca0012")
        #locations of things
        self.x_tail =  self.opti.variable(init_guess = 1, lower_bound = 0.7, upper_bound = 1.2)
        self.payload_mass_frac = self.opti.variable(init_guess = 10, lower_bound = 1, upper_bound = 5)
        # self.payload_mass_frac = 1

        #wing details
        self.masses = {}
        self.missions = ("M2", "M3")
        self.avl_mass = {mission: asb.MassProperties(mass=0) for mission in self.missions} #running totals for AVL mass file export, per mission
        self.total_mass = 0
        self.M2_max = 0.065
        self.M3_max = 28.16352056325056

        self.ratio_container_sensor = 0.5
        self.mass_sensor = self.opti.variable(init_guess = 1.5, lower_bound = 1, upper_bound = 3) #just the sensor
        self.mass_container = self.opti.variable(init_guess= 1, lower_bound = self.mass_sensor * self.ratio_container_sensor) #just the sensor

        #polar-based aero (optional): one polar per mission so each can use its own Reynolds number
        self.stall_margin = stall_margin
        self.flap_dCL = 1 #FLAPS: CL increment added to the polar's CL_max in the stall-speed check (takeoff/stall constraints)
        self.polar = None
        if polar_df is not None:
            Re_by_mission = Re if isinstance(Re, dict) else {m: Re for m in self.missions}
            self.polar = {m: AirfoilPolar(polar_df, Re_by_mission[m]) for m in self.missions}
        

    def _track_avl_mass(self, key, missions=None):
        #adds self.masses[key] into the running AVL mass total for each mission listed
        for mission in (missions if missions is not None else self.missions):
            self.avl_mass[mission] = self.avl_mass[mission] + self.masses[key]

    def optimize(self, verbose=True):   
        self.constants()        #establishes known constants
        self.run_dimensions()
        self.positions()
        self.run_aero_constraints()
        self.create_sensor()

        self.create_aero_plane()

        if self.given_airfoil:
            self.weights()          #sizes wing/spar/tail from the real airfoil geometry
        else:
            self.weights_no_af()    #generic 10% thick guess
        self.run_aero_no_af()       #uses the polar for CL/CD if one was given

        self.run_propulsion()
        self.run_structures()
        self.opti.subject_to(self.mass_empty > self.total_mass * 1.5)

        self.opti.maximize(self.objective())
        sol = self.opti.solve(verbose=verbose)
        if verbose:
            print(sol.value(self.mass_sensor))
        # for mission in self.missions:
        #   self.avl_mass[mission] = sol(self.avl_mass[mission]) #convert from symbolic opti variables to solved numeric values
        #   self.avl_mass[mission].export_AVL_mass_file(f"example_{mission}.mass")
        return sol

    def objective(self):
        mission = self.mission_no
        self.both_missions()
        if mission == 2:
            return self.M2()
        elif mission == 3:
            return self.M3()
        else:
            return self.both_missions()

#CONSTANT INITIALIZATION

    def constants(self):
        self.carbon_fiber_density = 1.75*1000 #kg/m^3
        self.carbon_fiber_yield_strength = 600e6 #Pa
        self.carbon_fiber_youngs_modulus = 230e9 #Pa
        self.carbon_fiber_shear_modulus = 50e9 #Pa
        self.carbon_fiber_layup_epoxy_factor = 2.2

        self.fiberglass_density = 2.6*1000 #kg/m^3
        self.fiberglass_areal_density = 0.017 #kg/m^2
        self.fiberglass_youngs_modulus = 70e9 #Pa
        self.fiberglass_shear_modulus = 30e9 #Pa
        self.fiberglass_layup_epoxy_factor = 2.2
        self.fiberglass_thickness = 0.0001

        self.kevlar_density = 1.4*1000 #kg/m^3
        self.kevlar_layup_epoxy_factor = 1 #per google, should fact check
        self.kevlar_youngs_modulus = 100e9 #Pa
        self.kevlar_yield_strength = 3e9 #Pa
        self.kevlar_shear_modulus = 717e6 #Pa

        self.basswood_density = 320 #kg/m^3
        self.basswood_youngs_modulus = 10e9 #Pa
        self.basswood_shear_modulus = 3.8e9 #Pa
        self.basswood_layup_epoxy_factor = 1.2
        self.basswoord_shear_strength = 6.8e6 #Pa

        self.plywood_density = 550 #kg/m^3
        self.plywood_thickness = 2.5/1000 #m

        self.balsa_density = 160 #kg/m^3
        self.balsa_youngs_modulus = 4e9 #Pa
        self.balsa_shear_modulus = 1.5e9 #Pa
        self.balsa_layup_epoxy_factor = 1.2
        self.balsa_shear_strength = 2.1e6 #Pa

        self.foam_density = 48 #kg/m^3

        self.monokote_density = 0.06 #kg/m^2
        self.carbon_fiber_areal_density = 0.2 #kg/m^2

        self.nu = 1.5111e-5 #m^2/s, kinematic viscosity of air at 68F
        self.rho = 1.10 #kg/m^3

        self.g = 9.81 #m/s^2

        self.lap_dist = 1000 #about 1km depending on how u do itLOL
        self.CL_max = 1.4
        self.S_g = 300 * units.foot

        self.N = 1 #cruise condition N = 1

#DIMENSION DERIVATION

    def run_dimensions(self):
        self.to_rad()
        self.wing_dimension()
        self.tail_dimension()
        
    def wing_dimension(self):   
        #b: span, c_r root chord, c_t tip chord, MAC mean aero chord
        self.b = np.sqrt(self.S * self.AR)
        self.mean_geo_chord = self.S / self.b
        self.c_r = 2 / (self.taper + 1) * self.mean_geo_chord
        self.c_t = self.taper * self.c_r
        self.MAC = 2/3 * self.c_r * (1 + self.taper + self.taper**2)/ (1 + self.taper)
        self.yMAC = self.b * (1 + 2 * self.taper) / (6 + 6 * self.taper)
        self.proj_span = self.b/np.cos(self.dihedral_rad)

    def tail_dimension(self):
        self.b_h = np.sqrt(self.S_h * self.AR_h)
        self.b_v = np.sqrt(self.S_v * self.AR_v)
        self.c_h = self.S_h / self.b_h
        self.c_v = self.S_v / self.b_v
        #rectangular tails

    def to_rad(self):
        self.dihedral_rad = self.dihedral * np.pi / 180
        self.washout_rad = self.washout * np.pi / 180
        self.decalage_rad = self.decalage * np.pi / 180

#AEROSANDBOX NATIVE AERODYNAMICS
    def create_wing(self):
        def const_quarter_chord(rc, tc):
            return rc / 4 - tc / 4

        tip_z = self.b / 2 * np.sin(self.dihedral * np.pi / 180)

        self.wing = asb.Wing(
            name="wing",
            symmetric=True,
            xsecs=[
                asb.WingXSec(
                    xyz_le=[0, 0, 0],  # UPDATED per Taras/AVL documentation
                    chord=self.c_r,
                    twist=0,
                    airfoil=self.airfoil,
                    control_surfaces=[
                        asb.ControlSurface(
                            name="Flap",
                            hinge_point=.6,
                            deflection=0,
                            symmetric=True
                        )
                    ]
                ),
                asb.WingXSec(
                    xyz_le=[
                        const_quarter_chord(self.c_r, self.c_t),
                        self.b / 4,
                        tip_z/2
                    ],  # UPDATED per Taras/AVL documentation
                    chord=self.c_r,
                    twist=0,
                    airfoil=self.airfoil,
                    control_surfaces=[
                        asb.ControlSurface(
                            name="Aileron",
                            hinge_point=.7,
                            deflection=0,
                            symmetric=False
                        )
                    ]
                ),
                asb.WingXSec(
                    xyz_le=[
                        const_quarter_chord(self.c_r, self.c_t),
                        self.b / 2,
                        tip_z
                    ],
                    chord=self.c_t,
                    twist=self.washout,
                    airfoil=self.airfoil
                )
            ]
        )


    def create_hor_stab(self):
        self.hor_stab = asb.Wing(
            name="horizontal_tail",
            symmetric=True,
            xsecs=[
                asb.WingXSec(
                    xyz_le=[0, 0, 0],
                    chord=self.c_h,
                    twist=self.decalage,
                    airfoil=self.airfoil_tail,
                    control_surfaces=[
                        asb.ControlSurface(
                            name="elevator",
                            hinge_point=0.70,
                            deflection=0
                        )
                    ]
                ),
                asb.WingXSec(
                    xyz_le=[0, self.b_h / 2, 0],
                    chord=self.c_h,
                    twist=self.decalage,
                    airfoil=self.airfoil_tail
                )
            ]
        ).translate([self.x_tail, 0, 0])

    def create_vert_stab(self):
        self.vert_stab = asb.Wing(
            name="vertical tail",
            xsecs=[
                asb.WingXSec(
                    xyz_le=[0, 0, 0],
                    chord=self.c_v,
                    twist=0,
                    airfoil=self.airfoil_tail,
                    control_surfaces=[
                        asb.ControlSurface(
                            name="rudder",
                            hinge_point=0.75,
                            deflection=0
                        )
                    ]
                ),
                asb.WingXSec(
                    xyz_le=[0, 0, self.b_v],
                    chord=self.c_v,
                    twist=0,
                    airfoil=self.airfoil_tail
                )
            ]
        ).translate([self.x_tail, 0, 0])

    def create_aero_plane(self):
        self.create_wing()
        self.create_hor_stab()
        self.create_vert_stab()
        self.plane = asb.Airplane(
            name="Prototype 0",
            xyz_ref=[self.c_r / 4, 0, 0],
            wings=[self.wing, self.hor_stab, self.vert_stab],
            fuselages = [self.fuselage, self.boom]
        )

    def run_aero_constraints(self):
        self.constraints_tail_dim()
        self.constraints_wing_dim()

    def constraints_wing_dim(self):
        self.opti.subject_to(self.b == 5.8 * units.foot)
        self.opti.subject_to(self.c_r < 0.5)

    def constraints_tail_dim(self):
        self.opti.subject_to(self.c_v == self.c_h)
        self.opti.subject_to(self.x_tail > self.c_r)

#ACTUAL OPERATING POINT ===============================================================================
    def op_point(self):
        self.velocity = {}
        self.alpha = {}
        self.q = {}
        self.op_pt = {}

        for mission in self.missions:
            self.velocity[mission] = self.opti.variable(
                init_guess=20,
                lower_bound=0,
                upper_bound=70 * units.mph
            )

            self.alpha[mission] = self.opti.variable(
                init_guess=0,
                upper_bound=13,
                lower_bound=-10
            )

            self.q[mission] = 1 / 2 * self.rho * self.velocity[mission]**2

            self.op_pt[mission] = asb.OperatingPoint(
                velocity=self.velocity[mission],
                alpha=self.alpha[mission],
                beta=0
            )

    
    def takeoff_constraint(self):
        if self.polar is not None:
            #this airfoil's own max lift instead of the fixed 1.4 / 2.0 guesses
            self.CL_max = self.polar["M2"].CL_max * self.stall_margin
        self.tw_to = (1.21 / (self.g * self.rho * self.CL_max * self.S_g) * self.mass["M2"] * self.g / self.S 
            + 0.05 ) #subject to mu and other stuff, look at 2025 design doc
        self.opti.subject_to(self.tw_to < 0.9)
        self.stall_speed = 15
        self.stall_CL = 1.5 if self.polar is None else self.polar["M2"].CL_max + self.flap_dCL
        self.opti.subject_to(1/2 * self.stall_speed ** 2 * self.rho * self.S * self.stall_CL > self.mass["M2"] * self.g)

    def stability_constraints(self):
        #static margin target and lateral/longitudinal stability derivatives, per mission
        #TODO: fuse_volume isn't modeled elsewhere yet - using the sensor pod's box volume as a stand-in.
        #swap this out if you build an actual fuselage shape, since a real body's volume (and its
        #moment-arm distribution) will differ from a plain box.
        self.fuse_volume = self.sensor_vol 

        self.static_margin = {}
        self.Cnb_corrected = {}
        self.Cma_corrected = {}
        for mission in self.missions:
            x_cg = self.avl_mass[mission].x_cg
            self.static_margin[mission] = (self.aero[mission]['x_np'] - x_cg) / self.wing.mean_aerodynamic_chord()
            self.Cnb_corrected[mission] = self.aero[mission]['Cnb'] - 2 * self.fuse_volume / (self.S * self.b) #corrected for fuselage contribution, which isn't in the basic plane model
            self.Cma_corrected[mission] = self.aero[mission]['Cma'] + 2 * self.fuse_volume / (self.S * self.wing.mean_aerodynamic_chord())

            self.opti.subject_to([
                self.static_margin[mission] == 0.15,
                self.Cnb_corrected[mission] > 0.04,
                self.Cnb_corrected[mission] < 0.1,
                self.Cma_corrected[mission] < -0.8,
                self.Cma_corrected[mission] > -1.8,
            ])

    def tail_volume(self):
        self.V_h = self.S_h *(self.x_tail + self.c_h/4)/ (self.S * self.c_r)
        self.V_v = self.S_v * (self.x_tail +self.c_h/4) / (self.S * self.b)
        self.opti.subject_to(self.V_h < 0.6)
        self.opti.subject_to(self.V_h > 0.3)
        self.opti.subject_to(self.V_v < 0.05)
        self.opti.subject_to(self.V_v > 0.02)
 

#REDUCED ORDER AERODYNAMICS 
    def run_aero_no_af(self):
        self.op_point()
        self.CL_no_af()
        self.CD_plane_no_af()
        self.drag_no_af()
        self.takeoff_constraint()
        self.tail_volume()

    def CL_no_af(self):
        #cruise
        # self.mass_empty = self.total_mass
        self.mass_empty = self.mass_sensor_total / self.payload_mass_frac


        #payload carried differs by mission: M2 carries sensor + container, M3 carries just the sensor
        self.payload_mass = {
            "M2": self.mass_sensor + self.mass_container,
            "M3": self.mass_sensor,
        }

        self.mass = {}
        
        for mission in self.missions:
            self.mass[mission] = self.mass_empty + self.payload_mass[mission]
            self.opti.subject_to(self.mass[mission] < 55 * units.pound)

            if self.polar is None:
                #placeholder: the CL that level flight needs, no airfoil involved
                self.CL[mission] = 2 * self.N * self.mass[mission] * self.g / (self.velocity[mission] ** 2 * self.S * self.rho)
            else:
                #CL comes from the polar at this mission's alpha; the optimizer solves alpha so that lift = weight
                polar = self.polar[mission]
                self.CL[mission] = polar.cl(self.alpha[mission])
                self.opti.subject_to([
                    self.CL[mission] * self.q[mission] * self.S / (self.N * self.mass[mission] * self.g) == 1,
                    self.alpha[mission] > polar.alpha_min,
                    self.alpha[mission] < polar.alpha_stall,
                    self.CL[mission] < self.stall_margin * polar.CL_max,
                ])

    def CD_plane_no_af(self):
        self.CDp = 0.002 #placeholder profile drag, only used when there's no polar
        self.e = 0.9 #oswald efficiency
        self.CD_factor = 2 #your original x2 on the wing CD - set to 1 to use raw (polar CD + induced)

        self.CDi = {}
        self.CD = {}
        self.CD_profile = {}
        for mission in self.missions:
            self.CDi[mission] = self.CL[mission] ** 2 / (self.e * np.pi * self.AR)
            if self.polar is None:
                self.CD_profile[mission] = self.CDp
            else:
                self.CD_profile[mission] = self.polar[mission].cd(self.alpha[mission])
            self.CD[mission] = self.CD_profile[mission] + self.CDi[mission]

    def drag_no_af(self):
        interference_factor = 1.08
        C_f = 0.004
        ratio = self.sensor_height / self.sensor_length  #fineness ratio (d/l) of the sensor pod as a streamlined body
        self.CD_fuse = 0.44 * ratio + 4 * C_f * (1 / ratio) + 4 * C_f * np.sqrt(ratio) #Hoerner 3-12, eqn 25 - on FRONTAL area, not S
        A_sensor_frontal = self.sensor_height ** 2

        A_fuse_frontal = 2 * self.sensor_height ** 2 
        self.CD_lg = 0.25
        A_lg = np.pi * 0.035**2 #frontal area of one gear leg

        #extra drag AREA (CD*A, not a CD-on-S) beyond the bare airframe, by mission - default 0, override as needed
        self.extra_drag_area = {mission: 0 for mission in self.missions}
        self.extra_drag_area["M3"] += self.CD_fuse * A_sensor_frontal #sensor pod exposed to the airstream on M3
        for mission in self.missions:
            self.extra_drag_area[mission] += self.CD_lg * A_lg #landing gear, assumed present on both missions
            self.extra_drag_area[mission] += self.CD_fuse * A_fuse_frontal

        self.drag = {}
        for mission in self.missions:
            self.drag[mission] = self.q[mission] * (self.CD[mission] * self.S + self.extra_drag_area[mission])
            self.drag[mission] *= 1.08

    def op_point_no_af(self):
        self.velocity = {}
        self.q = {}
        for mission in self.missions:
            self.velocity[mission] = self.opti.variable(
                init_guess=20,
                lower_bound=10,
                upper_bound=70 * units.mph
            )
            self.q[mission] = 1 / 2 * self.rho * self.velocity[mission]**2

#WEIGHT MODEL
        
    def positions(self):
        self.x_motor = self.opti.variable(init_guess = -.3, lower_bound= -0.4, upper_bound=-0.2)
        self.x_electronics = self.opti.variable(init_guess=-0.1, lower_bound=-0.5, upper_bound=4/5*self.c_r)
        self.x_battery = self.opti.variable(init_guess=-0.1, lower_bound=-0.5, upper_bound=4/5*self.c_r) #could update to vary between M2 and M3
        self.x_battery_avionics = self.opti.variable(init_guess=-0.2, lower_bound=-0.5, upper_bound=self.c_r)

    def misc_weight(self):
        self.mass_electronics = 0.130 + .0165 + .03 + .015 #130g esc (CC 100a) + 8ch elrs rx + 30 g of wire + 10a CC BEC (2-8s)
        self.mass_battery_avionics = 0.07 #2s 1Ah lipo
        self.mass_propeller = .56 #e.g. 16-18" prop and motor
        self.mass_tail_servo = 2*.035  #2x wing servos for tail
        self.mass_wing_servo = 4 * .035 #4x main wing servos
        self.mass_landing_gear = 0.2
        self.mass_misc = 1.0

        self.total_mass += self.mass_electronics + self.mass_battery_avionics + self.mass_propeller + self.mass_tail_servo + self.mass_wing_servo + self.mass_landing_gear + self.mass_misc

        self.masses['Electronics'] = asb.MassProperties(mass =self.mass_electronics , x_cg=self.x_electronics)
        self._track_avl_mass('Electronics')
        self.masses['Avionics battery'] = asb.MassProperties(mass=self.mass_battery_avionics, x_cg=self.x_battery_avionics) 
        self._track_avl_mass('Avionics battery')
        self.masses['Propeller'] = asb.MassProperties(mass=self.mass_propeller, x_cg = self.x_motor - 0.03) 
        self._track_avl_mass('Propeller')
        self.masses['Tail Servos'] = asb.MassProperties(mass=self.mass_tail_servo, x_cg = self.x_tail + self.c_h/2) 
        self._track_avl_mass('Tail Servos')
        self.masses['Wing servos'] = asb.MassProperties(mass=self.mass_wing_servo, x_cg = self.c_r/2) 
        self._track_avl_mass('Wing servos')
        self.masses['Landing Gear'] = asb.MassProperties(mass=self.mass_landing_gear, x_cg = 0) #UPDATED
        self._track_avl_mass('Landing Gear')
        self.masses['Misc'] = asb.MassProperties(mass=self.mass_misc, x_cg=0) #some of this will be the wing mount/interface, so probably around x=0
        self._track_avl_mass('Misc')

    def wing_weight(self): #not including the spar
        self.vol_foam = self.wing.volume() - units.inch * units.inch* self.b  # just the control surfaces, not that much right
        self.mass_foam = self.vol_foam * self.foam_density * self.fiberglass_layup_epoxy_factor
        self.mass_foam = self.mass_foam * 1.2

        self.mass_wing = self.mass_foam 
        self.total_mass += self.mass_wing 
        self.masses["Wing"] = asb.MassProperties(mass = self.mass_wing, x_cg = self.c_r/3) #check where c_g of wing would be
        self._track_avl_mass("Wing")

    def boom_weight(self):

        self.boom_weight = (self.x_tail - self.x_motor + self.c_h / 4) / 1.27 * 0.173 #scaling from last year
        self.total_mass += self.boom_weight
        self.masses["Boom"] = asb.MassProperties(mass = self.boom_weight, x_cg = (self.x_tail + self.c_h/4 + self.x_motor) / 2 )
        self._track_avl_mass("Boom")

    def spar_weight(self):
        
        self.spar_cap_thickness = 0.002 #2mm worth of carbon?
        self.spar_wood_thickness = units.inch
        self.spar_wood_width = units.inch
        self.vol_spar_cap = 4 * self.spar_cap_thickness * self.proj_span * self.spar_wood_width 
        self.mass_spar_cap = self.vol_spar_cap * self.carbon_fiber_density * self.fiberglass_layup_epoxy_factor
        self.mass_spar_misc = 0.01
        self.mass_spar = self.mass_spar_cap + self.mass_spar_misc

        self.total_mass += self.mass_spar
        self.masses["Spar"] = asb.MassProperties(mass = self.mass_spar, x_cg = self.c_r / 4)
        self._track_avl_mass("Spar")

    def spar_weight_no_af(self):
        
        self.spar_cap_thickness = 0.002 #2mm worth of carbon?
        self.spar_wood_thickness = units.inch
        self.spar_wood_width = units.inch
        self.vol_spar_cap = 4 * self.spar_cap_thickness * self.proj_span * self.spar_wood_width 
        self.mass_spar_cap = self.vol_spar_cap * self.carbon_fiber_density * self.fiberglass_layup_epoxy_factor
        self.mass_spar_misc = 0.01
        self.mass_spar = self.mass_spar_cap + self.mass_spar_misc

        self.total_mass += self.mass_spar
        self.masses["Spar"] = asb.MassProperties(mass = self.mass_spar, x_cg = self.c_r / 4)
        self._track_avl_mass("Spar")

    def tail_weight(self):
        self.vol_tail = self.hor_stab.volume() + self.vert_stab.volume()
        self.mass_foam_tail = self.vol_tail * self.foam_density

        self.area_wetted_tail = self.hor_stab.area(type="wetted") + self.vert_stab.area(type="wetted")
        self.mass_wetted_tail = self.area_wetted_tail * self.fiberglass_areal_density * self.fiberglass_layup_epoxy_factor

        self.mass_tail = self.mass_foam_tail + self.mass_wetted_tail

        self.total_mass += self.mass_tail
        self.masses["Tail"] = asb.MassProperties(mass = self.mass_tail, x_cg = self.x_tail + self.c_h / 3)
        self._track_avl_mass("Tail")

    def wing_weight_no_af(self): #not including the spar
        self.area_chord = 0.064 #from xfoil
        
        self.n_stringer = 10
        self.stringer_width = 0.007
        self.stringer_thickness = 0.001
        self.vol_stringer_wood = self.proj_span * self.n_stringer * self.stringer_width * self.stringer_thickness 

        self.n_ribs = 7 #for one half of the plane
        self.rib_chords = 2 * [self.c_r + (self.c_t - self.c_r) * i / (self.n_ribs - 1) for i in range(self.n_ribs)]
        self.vol_wood = (0.75 * self.plywood_thickness * self.area_chord * sum(c**2 for c in self.rib_chords))
        self.mass_plywood = self.plywood_density * self.vol_wood
        self.mass_stringer = self.vol_stringer_wood * self.basswood_density
        self.mass_wood = self.mass_plywood + self.mass_stringer
        self.mass_wood = self.mass_wood * 1.35 #correction factor for extra glue

        self.vol_foam = 0.2 * self.area_chord * self.b * (self.c_r)**2 # just the control surfaces, not that much right
        self.vol_control_wood = 0.08 * self.c_r * self.plywood_thickness * self.proj_span
        self.mass_control_wood = self.vol_control_wood * self.plywood_density
        self.mass_foam = self.vol_foam * self.foam_density + self.mass_control_wood
        self.mass_foam = self.mass_foam * 1.2

        self.skin_monokote = self.c_r * (0.7 - 0.25) * 2 * 1.1 * self.proj_span * self.monokote_density #1.1 is correction for curvature
        self.skin_carbon_fiber = self.c_r * 0.25 * 1.3 * 2 * self.proj_span * self.carbon_fiber_areal_density #1.3 is correction for curvature
        self.mass_skin = self.skin_monokote + self.skin_carbon_fiber
        self.mass_skin = self.mass_skin * 1.1 #1.1 correction factor

        self.mass_wing = self.mass_wood + self.mass_foam + self.mass_skin

        self.total_mass += self.mass_wing
        self.masses["Wing"] = asb.MassProperties(mass = self.mass_wing, x_cg = self.c_r/3) #check where c_g of wing would be
        self._track_avl_mass("Wing")

    def tail_weight_no_af(self):
        self.area_chord_tail = 0.07 #a guess for area of chord length 1

        self.vol_tail = self.area_chord_tail * self.c_h ** 2 * (self.b_h + self.b_v)
        self.mass_foam_tail = self.vol_tail * self.foam_density

        self.area_wetted_tail = self.c_v * 2.2 * self.b_v + self.c_h * 2.2 * self.b_h #this is my guess for area two sides, around chord length 1.1x per side?
        self.mass_wetted_tail = self.area_wetted_tail * self.fiberglass_areal_density * self.fiberglass_layup_epoxy_factor

        self.mass_tail = self.mass_foam_tail + self.mass_wetted_tail
        self.total_mass += self.mass_tail
        self.masses["Tail"] = asb.MassProperties(mass = self.mass_tail, x_cg = self.x_tail + self.c_h / 3)
        self._track_avl_mass("Tail")

    def weights(self):
        self.wing_weight()
        self.misc_weight()
        self.tail_weight()
        self.boom_weight()
        self.spar_weight()

    def weights_no_af(self):
        self.positions()
        self.wing_weight_no_af()
        self.misc_weight()
        self.tail_weight_no_af()
        self.boom_weight()
        self.spar_weight_no_af()

#MISSION DEFINITIONS

    def M2(self):
        self.M2_score = self.mass_sensor_total / self.duration["M2"]
        self.opti.subject_to(self.duration["M2"]<300)
        return self.M2_score    

    def M3(self): 
        self.M3_laps = self.velocity["M3"] * 5 * 60 / self.lap_dist
        self.M3_score = self.mass_sensor * self.M3_laps
        return self.M3_score

    def both_missions(self):
        self.both_mission_score = self.M2() / self.M2_max + self.M3() / self.M3_max
        return self.both_mission_score

#PROPULSION VAGUE

    def run_propulsion(self):
        self.general_power()

        #flight duration differs by mission: M2 is however long 5 laps s at its solved speed,
        #M3 is a fixed 5 minute window - override per mission as needed
        self.M2_distance = self.lap_dist * 5 # 5 laps
        self.extra_battery = {}
        self.duration = {
            "M2": self.safety_factor * self.M2_distance / self.velocity["M2"],
            "M3": self.safety_factor * 5 * 60,
        }
        
        self.power = self.drag["M2"] * self.velocity["M2"] / 0.6
        self.opti.subject_to(self.max_amperage * self.voltage > self.power)

        for mission in self.missions:
            self.opti.subject_to(self.drag[mission] * self.velocity[mission] * self.duration[mission] < self.battery_power )
            self.extra_battery[mission] = (self.battery_power- self.drag[mission] * self.velocity[mission] * self.duration[mission]) / self.battery_power  
        
    def general_power(self):
        self.battery_power = 100 * 3600 #J
        self.voltage = 25
        self.max_amperage = 100
        self.safety_factor = 1.8

#SENSOR INITIALIZATION

    def sensor_dim(self):
        self.sensor_height = self.opti.variable(init_guess = 0.12, lower_bound = 5 * units.inch, upper_bound = 10 * units.inch)
        self.sensor_width = self.opti.variable(init_guess = 0.1, lower_bound = 3 * units.inch, upper_bound = 6 * units.inch)
        self.sensor_length = self.opti.variable(init_guess = 12 * units.inch, lower_bound = 10 * units.inch, upper_bound = 14 * units.inch)

        self.sensor_vol = self.sensor_length * self.sensor_width * self.sensor_height


        self.fuselage = asb.Fuselage(
            xsecs=[
                asb.FuselageXSec(
                    xyz_c=[-0.06, 0,  -self.sensor_height/1.1],
                    radius=0,
                ),
                asb.FuselageXSec(
                    xyz_c=[-0.04, 0,  -self.sensor_height/1.1],
                    radius=self.sensor_height/2,
                ),
                asb.FuselageXSec(
                    xyz_c=[0, 0, -self.sensor_height/1.1],
                    radius=self.sensor_height,
                ),
                asb.FuselageXSec(
                    xyz_c=[self.sensor_length*1.5, 0, -self.sensor_height/1.1],
                    radius=self.sensor_height,
                ),
                asb.FuselageXSec(
                    xyz_c=[self.sensor_length*1.5 + 0.04, 0,  -self.sensor_height/1.1],
                    radius=self.sensor_height/2,
                ),
                asb.FuselageXSec(
                    xyz_c=[self.sensor_length*1.5 + 0.06, 0,  -self.sensor_height/1.1],
                    radius=0,
                ),

            ]
        )

        self.boom = asb.Fuselage(
            xsecs=[
                asb.FuselageXSec(
                    xyz_c=[self.x_motor, 0,  0],
                    radius=0,
                ),
                asb.FuselageXSec(
                    xyz_c=[self.x_motor+ 0.02, 0, 0],
                    radius=0.0254,
                ),
                asb.FuselageXSec(
                    xyz_c=[self.x_tail,  0,0],
                    radius=0.0254,
                ),
                asb.FuselageXSec(
                    xyz_c=[self.x_tail + 0.02, 0,0],
                    radius=0,
                ),

            ]
        )

    def sensor_weight(self):
        self.mass_sensor_total = self.mass_sensor + self.mass_container
        # self.opti.subject_to(self.mass_sensor_total < 12)

        self.x_sensor = self.opti.variable(init_guess=0, lower_bound = -0.1, upper_bound=0.2) #TODO: set this to the actual sensor/container CG location, this is just a placeholder
        self.masses['Sensor'] = asb.MassProperties(mass = self.mass_sensor, x_cg = self.x_sensor)
        self.masses['Container'] = asb.MassProperties(mass = self.mass_container, x_cg = self.x_sensor)

        #M2 flies with both the sensor and its container onboard; M3 flies with just the sensor
        self._track_avl_mass('Sensor') #both missions
        self._track_avl_mass('Container', missions=('M2',)) #M2 only

    def sensor_constraint(self):
        self.opti.subject_to((self.mass_sensor + self.mass_container) / (self.sensor_vol) < 3000) #less than straight steel

    def create_sensor(self):
        self.sensor_dim()
        self.sensor_weight()
        self.sensor_constraint()

#STRUCTURE INITIALIZATION

    def run_structures(self):
        self.wing_bending()

    def torsion(opti, params):
        #checks the torsion of the tail to the structures of the wing
        pass
        
    def wing_bending(self):
        self.spar_I = 2 * (self.spar_wood_thickness * self.spar_cap_thickness**3 / 12 + (self.spar_wood_thickness/2) **2 * self.spar_cap_thickness * self.spar_wood_width)
        if not self.given_airfoil:
            self.spar_I /= 5
        self.N_max = 7
        self.max_load_per_length_M2 = self.N_max * self.g * self.mass["M2"] / self.proj_span

        self.deflection_M2 = self.max_load_per_length_M2 * (self.b / 2)**4 / (8 * self.carbon_fiber_youngs_modulus * self.spar_I)
        self.deflection_max = 0.05
        self.opti.subject_to(self.deflection_M2 < self.deflection_max)

#print statements
    

    def get_results(self, sol):
        #pulls the solved values into a plain nested dict of real numbers - safe to print or save as JSON
        return {
            "geometry": {
                "wing_area_m2": sol(self.S),
                "aspect_ratio": sol(self.AR),
                "span_m": sol(self.b),
                "root_chord_m": sol(self.c_r),
                "mean_aero_chord_m": sol(self.MAC),
            },
            "mass_kg": {
                "empty": sol(self.mass_empty),
                "wing": sol(self.mass_wing),
                "tail": sol(self.mass_tail),
                "sensor": sol(self.mass_sensor),
                "container": sol(self.mass_container),
                "total_M2": sol(self.mass["M2"]),
                "total_M3": sol(self.mass["M3"]),
            },
            "missions": {
                mission: {
                    "velocity_mps": sol(self.velocity[mission]),
                    "CL": sol(self.CL[mission]),
                    # "CD": sol(self.CD[mission]),
                    "drag_N": sol(self.drag[mission]),
                }
                for mission in self.missions
            },
            "scores": {
                "M2_score": sol(self.M2_score),
                "M3_score": sol(self.M3_score),
            },
        }

    def print_summary(self, sol):
        results = self.get_results(sol)
        print("\n" + "=" * 50)
        print("OPTIMIZATION RESULT".center(50))
        print("=" * 50)
        for section, values in results.items():
            print(f"\n-- {section} --")
            if section == "missions":
                for mission, vals in values.items():
                    print(f"  {mission}:")
                    for k, v in vals.items():
                        print(f"    {k:<18} {v:>10.4f}")
            else:
                for k, v in values.items():
                    print(f"  {k:<20} {v:>10.4f}")
        print("=" * 50 + "\n")

    def save_solution(self, sol, path="solution.json"):
        results = self.get_results(sol)
        with open(path, "w") as f:
            json.dump(results, f, indent=2)
        print(f"Saved solution to {path}")


# =====================================================================================================
# SWEEP: run the optimizer once per airfoil
# =====================================================================================================

def run_airfoil(name, polar_df, Re=500_000, mission=0, verbose=False):
    """
    Build and solve the optimizer for ONE airfoil. Always returns a flat dict (never raises),
    so one bad airfoil can't kill a long sweep.  status is "ok" or "failed" (see the `error` column).
    """
    row = {"airfoil": name, "status": "ok", "error": ""}
    t0 = time.time()
    try:
        v = variables(airfoil=name, polar_df=polar_df, Re=Re)
        v.mission_no = mission
        sol = v.optimize(verbose=verbose)
        res = v.get_results(sol)

        score = {2: v.M2_score, 3: v.M3_score}.get(mission, v.both_mission_score)
        row["score"] = float(sol(score))
        for section in ("geometry", "mass_kg", "scores"):
            row.update({k: float(x) for k, x in res[section].items()})
        for m in v.missions:
            row.update({f"{k}_{m}": float(x) for k, x in res["missions"][m].items()})
            row[f"alpha_deg_{m}"] = float(sol(v.alpha[m]))
            row[f"CD_profile_{m}"] = float(sol(v.CD_profile[m]))
            row[f"L_over_D_{m}"] = row[f"CL_{m}"] / float(sol(v.CD[m]))
            #Reynolds number the solution actually flies at vs the one the polar was read at
            row[f"Re_assumed_{m}"] = v.polar[m].Re
            row[f"Re_actual_{m}"] = float(sol(v.velocity[m]) * sol(v.MAC) / v.nu)
        row["t_over_c"] = float(v.airfoil.max_thickness())
        row["CL_max_polar"] = v.polar["M2"].CL_max
        row["Re_clamped"] = any(p.Re_clamped for p in v.polar.values()) #True = the data didn't cover the Re you asked for
    except Exception as e:
        msg = str(e)
        status = re.search(r"return_status is '(\w+)'", msg)
        row["status"] = "failed"
        row["error"] = status.group(1) if status else msg.strip().splitlines()[-1][:150]
    row["time_s"] = round(time.time() - t0, 2)
    return row


def _run_job(job):  # top-level so multiprocessing can pickle it
    return run_airfoil(*job)

# =====================================================================================================
# SWEEP: find maximum airfoil L/D at a specified Reynolds number
# =====================================================================================================

def max_ld_airfoil(name, polar_df, Re=500_000):
    """
    Find the maximum 2D airfoil L/D for one airfoil at the requested Reynolds number.

    The polar is interpolated in log(Re) using AirfoilPolar, and only the
    pre-stall branch is considered.
    """
    row = {
        "airfoil": name,
        "Re_requested": Re,
        "status": "ok",
        "error": "",
    }

    try:
        polar = AirfoilPolar(polar_df, Re)

        # Values from the interpolated pre-stall polar
        alpha = polar.alpha_grid
        CL = polar.CL_grid
        CD = polar.CD_grid

        # Airfoil L/D
        LD = CL / (CD + CL ** 2 / (2 * np.pi * 4 * 0.8))

        # Maximum L/D
        i_max = onp.nanargmax(LD)

        row["Re_actual"] = polar.Re
        row["Re_clamped"] = polar.Re_clamped

        row["alpha_deg"] = alpha[i_max]
        row["CL"] = CL[i_max]
        row["CD"] = CD[i_max]
        row["L_over_D"] = LD[i_max]

        row["CL_max"] = polar.CL_max
        row["alpha_stall_deg"] = polar.alpha_stall
        row["t_over_c"] = asb.Airfoil(name.lower()).max_thickness()

    except Exception as e:
        row["status"] = "failed"
        row["error"] = str(e)[:200]

    return row


def sweep_airfoils(airfoils, Re=500_000, out_csv="max_LD_results.csv"):
    """
    Find maximum 2D airfoil L/D for every airfoil at one Reynolds number.

    Saves results to CSV and returns a DataFrame sorted by maximum L/D.
    """

    rows = []

    for i, (name, polar_df) in enumerate(airfoils.items(), 1):

        print(f"[{i}/{len(airfoils)}] {name}")

        row = max_ld_airfoil(
            name=name,
            polar_df=polar_df,
            Re=Re,
        )

        rows.append(row)

    df = pd.DataFrame(rows)

    # Put successful airfoils first and sort by L/D
    ok = df[df["status"] == "ok"].sort_values(
        "L_over_D",
        ascending=False
    )

    failed = df[df["status"] != "ok"]

    df = pd.concat([ok, failed], ignore_index=True)

    df.to_csv(out_csv, index=False)

    print(f"\nSaved results to {out_csv}")

    print("\nTop 15 airfoils:")
    print(
        df[df.status == "ok"]
        .head(15)
        [
            [
                "airfoil",
                "Re_actual",
                "alpha_deg",
                "CL",
                "CD",
                "L_over_D",
                "CL_max",
            ]
        ]
        .to_string(index=False)
    )

    if len(failed):
        print(f"\n{len(failed)} airfoils failed:")
        print(df[df.status != "ok"][["airfoil", "error"]].to_string(index=False))

    return df

def sweep(airfoils, Re=500_000, mission=0, names=None, n_jobs=1, out_csv="sweep_results.csv", resume=True):
    """
    Run the optimizer for every airfoil in `airfoils` (the dict from load_airfoils).
    Re       : Reynolds number for ALL airfoils (or {"M2": .., "M3": ..})
    mission  : 2 = M2 only, 3 = M3 only, anything else = combined score
    n_jobs   : processes to run in parallel (each solve is single-threaded, so this scales well)
    resume   : if out_csv exists, skip airfoils already in it (so you can stop/restart a long run)
    Returns a DataFrame, best score first (failed airfoils at the bottom).
    """
    from concurrent.futures import ProcessPoolExecutor
    names = list(names) if names is not None else list(airfoils)

    rows = []
    if resume and Path(out_csv).exists():
        prev = pd.read_csv(out_csv)
        rows = prev.to_dict("records")
        names = [n for n in names if n not in set(prev["airfoil"])]
        print(f"Resuming: {len(prev)} airfoils already in {out_csv}, {len(names)} left")

    jobs = [(n, airfoils[n], Re, mission) for n in names]
    t0 = time.time()

    def checkpoint(i):
        pd.DataFrame(rows).to_csv(out_csv, index=False)
        n_ok = sum(r["status"] == "ok" for r in rows)
        print(f"  {i}/{len(jobs)} done  ({n_ok} ok, {len(rows) - n_ok} failed)  {time.time() - t0:.0f}s elapsed")

    if n_jobs == 1:
        for i, job in enumerate(jobs, 1):
            rows.append(run_airfoil(*job))
            if i % 25 == 0:
                checkpoint(i)
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            for i, row in enumerate(pool.map(_run_job, jobs, chunksize=4), 1):
                rows.append(row)
                if i % 25 == 0:
                    checkpoint(i)
    checkpoint(len(jobs))

    df = pd.DataFrame(rows)
    ok = df[df.status == "ok"].sort_values("score", ascending=False)
    return pd.concat([ok, df[df.status != "ok"]], ignore_index=True)

def main2():
    import argparse

    p = argparse.ArgumentParser()

    p.add_argument(
        "--folder",
        required=True,
        help="folder with the <airfoil>_Re<number>.txt polar files"
    )

    p.add_argument(
        "--Re",
        type=float,
        default=1_000_000,
        help="Reynolds number at which to evaluate every airfoil"
    )

    p.add_argument(
        "--out",
        default="max_LD_results.csv",
        help="output CSV filename"
    )

    args = p.parse_args()

    # Load all airfoil polars
    airfoils = load_airfoils(args.folder)

    print(f"Loaded {len(airfoils)} airfoils")
    print(f"Evaluating maximum L/D at Re = {args.Re:,.0f}\n")

    results = sweep_airfoils(
        airfoils,
        Re=args.Re,
        out_csv=args.out,
    )

def main1():
    import argparse
    import casadi as ca

    try:
        ca.GlobalOptions.setNumpyMode(-1)
    except AttributeError:
        pass

    p = argparse.ArgumentParser()
    p.add_argument("--folder", required=True, help="folder with the <airfoil>_Re<number>.txt polars")
    p.add_argument("--Re", type=float, default=500000, help="Reynolds number used for every airfoil")
    p.add_argument("--mission", type=int, default=0, help="2 = M2 only, 3 = M3 only, anything else = combined")
    p.add_argument("--airfoil", default=None, help="run just this airfoil (full solver output) instead of sweeping all")
    p.add_argument("--jobs", type=int, default=1, help="parallel processes for the sweep")
    p.add_argument("--out", default="sweep_results.csv")
    args = p.parse_args()

    airfoils = load_airfoils(args.folder)
    print(f"Loaded {len(airfoils)} airfoils")

    if args.airfoil:
        v = variables(airfoil=args.airfoil, polar_df=airfoils[args.airfoil], Re=args.Re)
        v.mission_no = args.mission
        sol = v.optimize(verbose=True)
        v.print_summary(sol)
        v.save_solution(sol, f"solution_{args.airfoil}.json")
    else:
        results = sweep(airfoils, Re=args.Re, mission=args.mission, n_jobs=args.jobs, out_csv=args.out)
        show = ["airfoil", "score", "M2_score", "M3_score", "wing_area_m2", "aspect_ratio", "t_over_c", "CL_max_polar", "time_s"]
        print("\nTop 15:")
        print(results[results.status == "ok"].head(15)[show].to_string(index=False))
        failed = results[results.status != "ok"]
        if len(failed):
            print(f"\n{len(failed)} failed - reasons:")
            print(failed.error.value_counts().to_string())

if __name__ == "__main__":
    main1()