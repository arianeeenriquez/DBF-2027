import aerosandbox as asb
import aerosandbox.numpy as np
from aerosandbox.tools import units


class variables:
	def __init__(self):
		self.opti = asb.Opti()
		self.order = 2

		#overall wing design
		self.S = self.opti.variable(init_guess = 0.1, lower_bound = 0.01, upper_bound = 1)
		self.AR = self.opti.variable(init_guess = 6, lower_bound = 4, upper_bound = 15)
		self.taper = 1
		self.airfoil = asb.Airfoil("sd7032")

		self.dihedral = 5 #degrees
		self.washout = 0 #degrees
		#tail design
		self.AR_v = self.opti.variable(init_guess = 2, lower_bound = 1, upper_bound = 5)
		self.AR_h = self.opti.variable(init_guess = 4, lower_bound = 1, upper_bound = 5)
		self.S_v = self.opti.variable(init_guess = 0.01, lower_bound = 0.001, upper_bound = .1)
		self.S_h = self.opti.variable(init_guess = 0.01, lower_bound = 0.001, upper_bound = .1)

		self.decalage = 5 #degree

		self.airfoil_tail = asb.Airfoil("naca0012")
		#locations of things
		self.x_tail =  self.opti.variable(init_guess = 1, lower_bound = 0.5, upper_bound = 1.2)

		#wing details
		self.masses = {}
		self.total_mass = 0
		self.M2_max = 0.48
		self.M3_max = 2.4

	def optimize(self):	
		constants(self) 		#establishes known constants
		run_dimensions(self)
		run_aero_constraints(self)
		create_sensor(self)

		given_airfoil = False
		if given_airfoil:
			#creates a aerosandbox entity of the wing, hor stab and vert stab
			create_wing(self)
			create_hor_stab(self)
			create_vert_stab(self)

			#greates a mass model based on the wing and tail dimensions
			weights(self)
		else:
			weights_no_af(self)

		run_aero(self)
		run_propulsion(self)
		run_structures(self)

		self.opti.maximize(self.objective())
		sol = self.opti.solve()
		return sol

	def objective(self):
		return both_missions(self)

#CONSTANT INITIALIZATION

	def constants(v):
		v.carbon_fiber_density = 1.75*1000 #kg/m^3
		v.carbon_fiber_yield_strength = 600e6 #Pa
		v.carbon_fiber_youngs_modulus = 230e9 #Pa
		v.carbon_fiber_shear_modulus = 50e9 #Pa
		v.carbon_fiber_layup_epoxy_factor = 2.2

		v.fiberglass_density = 2.6*1000 #kg/m^3
		v.fiberglass_areal_density = 0.017 #kg/m^2
		v.fiberglass_youngs_modulus = 70e9 #Pa
		v.fiberglass_shear_modulus = 30e9 #Pa
		v.fiberglass_layup_epoxy_factor = 2.2
		v.fiberglass_thickness = 0.0001

		v.kevlar_density = 1.4*1000 #kg/m^3
		v.kevlar_layup_epoxy_factor = 1 #per google, should fact check
		v.kevlar_youngs_modulus = 100e9 #Pa
		v.kevlar_yield_strength = 3e9 #Pa
		v.kevlar_shear_modulus = 717e6 #Pa

		v.basswood_density = 320 #kg/m^3
		v.basswood_youngs_modulus = 10e9 #Pa
		v.basswood_shear_modulus = 3.8e9 #Pa
		v.basswood_layup_epoxy_factor = 1.2
		v.basswoord_shear_strength = 6.8e6 #Pa

		v.plywood_density = 550 #kg/m^3
		v.plywood_thickness = 2.5/1000 #m

		v.balsa_density = 160 #kg/m^3
		v.balsa_youngs_modulus = 4e9 #Pa
		v.balsa_shear_modulus = 1.5e9 #Pa
		v.balsa_layup_epoxy_factor = 1.2
		v.balsa_shear_strength = 2.1e6 #Pa

		v.foam_density = 48 #kg/m^3

		v.monokote_density = 0.06 #kg/m^2
		v.carbon_fiber_areal_density = 0.2 #kg/m^2

		v.nu = 1.5111e-5 #m^2/s, kinematic viscosity of air at 68F
		v.rho = 1.10 #kg/m^3

		v.g = 9.81 #m/s^2

		v.lap_dist = 1000 #about 1km depending on how u do itLOL

#DIMENSION DERIVATION

	def run_dimensions(self):
		to_rad(self)
		wing_dimension(self)
		tail_dimension(self)
		
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

#REDUCED ORDER AERODYNAMICS	
	def run_aero_no_af(v):
		op_point(v)
		CL_no_af(v)
		CD_plane_no_af(v)
		drag_no_af(v)

	def CL_no_af(v):
		#cruise
		v.mass_empty = v.total_mass
		v.mass_M2 = v.mass_empty + v.mass_sensor + v.mass_container
		v.mass_M3 = v.mass_empty = v.mass_sensor

		v.opti.subject_to(v.mass_M2 < 55 * units.pound)
		v.opti.subject_to(v.mass_M3 < 55 * units.pound)

		v.N = 1 #cruise condition N = 1
		v.CL_M2 = 2 * v.N * v.mass_M2 * v.g /(v.velocity_M2 ** 2 * v.S * v.rho) 
		v.CL_M3 = 2 * v.N * v.mass_M3 * v.g /(v.velocity_M3 ** 2 * v.S * v.rho) 

		v.opti.subject_to(v.CL_M2 < 0.5) #cruise CL
		v.opti.subject_to(v.CL_M3 < 0.5)

	def CD_plane_no_af(v):
		v.CDp = 0.030 
		v.e = 0.9 #oswald efficiency
		v.CDi_M2 = v.CL_M2 ** 2 / (v.e * np.pi * v.AR)
		v.CD_M2 = v.CDp + v.CDi_M2

		v.CDi_M3 = v.CL_M3 ** 2 / (v.e * np.pi * v.AR) 
		v.CD_M3 = v.CDp + v.CDi_M2

	def drag_no_af(v):
		v.CD_sensor = 0.02

		v.drag_M2 = v.CD_M2 * v.q_M2 * v.S
		v.drag_M3 = v.q_M3 * (v.CD_M3 * v.S + v.CD_sensor * v.sensor_width * v.sensor_height)

	def op_point_no_af(v):
	    v.velocity_M2 = v.opti.variable(
	        init_guess=20,
	        lower_bound=10,
	        upper_bound=85 * units.mph
	    )

	    v.velocity_M3 = v.opti.variable(
	        init_guess=20,
	        lower_bound=10,
	        upper_bound=85 * units.mph
	    )

	    v.q_M2 = 1 / 2 * v.rho * v.velocity_M2**2
	    v.q_M3 = 1 / 2 * v.rho * v.velocity_M3**2

#WEIGHT MODEL
		
	def positions(v):
		v.x_motor = v.opti.variable(init_guess = -.3, lower_bound= -0.5, upper_bound=0)
		v.x_electronics = v.opti.variable(init_guess=-0.1, lower_bound=-0.5, upper_bound=4/5*v.c_r)
		v.x_battery = v.opti.variable(init_guess=-0.1, lower_bound=-0.5, upper_bound=4/5*v.c_r) #could update to vary between M2 and M3
		v.x_battery_avionics = v.opti.variable(init_guess=-0.2, lower_bound=-0.5, upper_bound=v.c_r)

	def misc_weight(v):
		v.mass_electronics = 0.130 + .0165 + .03 + .015 #130g esc (CC 100a) + 8ch elrs rx + 30 g of wire + 10a CC BEC (2-8s)
		v.mass_battery_avionics = 0.07 #2s 1Ah lipo
		v.mass_propeller = .06 #e.g. 16-18" prop
		v.mass_tail_servo = 2*.035  #2x wing servos for tail
		v.mass_wing_servo = 4 * .035 #4x main wing servos
		v.mass_landing_gear = 0.2
		v.mass_misc = 0.4

		v.total_mass += v.mass_electronics + v.mass_battery_avionics + v.mass_propeller + v.mass_tail_servo + v.mass_wing_servo + v.mass_landing_gear + v.mass_misc

		v.masses['Electronics'] = asb.MassProperties(mass =v.mass_electronics , x_cg=v.x_electronics)
		v.masses['Avionics battery'] = asb.MassProperties(mass=v.mass_battery_avionics, x_cg=v.x_battery_avionics) 
		v.masses['Propeller'] = asb.MassProperties(mass=v.mass_propeller, x_cg = v.x_motor - 0.03) 
		v.masses['Tail Servos'] = asb.MassProperties(mass=v.mass_tail_servo, x_cg = v.x_tail + v.c_h/2) 
		v.masses['Wing servos'] = asb.MassProperties(mass=v.mass_wing_servo, x_cg = v.c_r/2) 
		v.masses['Landing Gear'] = asb.MassProperties(mass=v.mass_landing_gear, x_cg = 0) #UPDATED
		v.masses['Misc'] = asb.MassProperties(mass=v.mass_misc, x_cg=0) #some of this will be the wing mount/interface, so probably around x=0

	def wing_weight(v): #not including the spar
		v.area_chord = 6.43354654E-02 #from xfoil
		
		v.n_stringer = 10
		v.stringer_width = 0.007
		v.stringer_thickness = 0.001
		v.vol_stringer_wood = v.proj_span * v.n_stringer * v.stringer_width * v.stringer_thickness 

		v.n_ribs = 7 #for one half of the plane
		v.rib_chords = 2 * [v.c_r + (v.c_t - v.c_r) * i / (v.n_ribs - 1) for i in range(v.n_ribs)]
		v.vol_wood = (0.75 * v.plywood_thickness * v.area_chord * sum(c**2 for c in v.rib_chords))
		v.mass_plywood = v.plywood_density * v.vol_wood
		v.mass_stringer = v.vol_stringer_wood * v.basswood_density
		v.mass_wood = v.mass_plywood + v.mass_stringer
		v.mass_wood = v.mass_wood * 1.35 #correction factor for extra glue

		v.vol_foam = 0.2 * v.wing.volume()  # just the control surfaces, not that much right
		v.vol_control_wood = v.airfoil.local_thickness(x_over_c = 0.6) * v.c_r * v.plywood_thickness * v.proj_span
		v.mass_control_wood = v.vol_control_wood * v.plywood_density
		v.mass_foam = v.vol_foam * v.foam_density + v.mass_control_wood
		v.mass_foam = v.mass_foam * 1.2

		v.skin_monokote = v.c_r * (0.7 - 0.25) * 2 * 1.1 * v.proj_span * v.monokote_density #1.1 is correction for curvature
		v.skin_carbon_fiber = v.c_r * 0.25 * 1.3 * 2 * v.proj_span * v.carbon_fiber_areal_density #1.3 is correction for curvature
		v.mass_skin = v.skin_monokote + v.skin_carbon_fiber
		v.mass_skin = v.mass_skin * 1.1 #1.1 correction factor

		v.mass_wing = v.mass_wood + v.mass_foam + v.mass_skin
		v.total_mass += v.mass_wing
		v.masses["Wing"] = asb.MassProperties(mass = v.mass_wing, x_cg = v.c_r/3) #check where c_g of wing would be

	def boom_weight(v):
		v.boom_weight = (v.x_tail - v.x_motor + v.c_h / 4) / 1.27 * 0.173 #scaling from last year
		v.total_mass += v.boom_weight
		v.masses["Boom"] = asb.MassProperties(mass = v.boom_weight, x_cg = (v.x_tail + v.c_h/4 + v.x_motor) / 2 )

	def spar_weight(v):
		
		v.spar_cap_thickness = 0.0004 #2mm worth of carbon?
		v.spar_wood_thickness = v.airfoil.local_thickness(x_over_c=0.25) * v.c_r - v.spar_cap_thickness
		v.spar_wood_width = 0.0127 #1/2 inch balsa for now
		v.vol_spar_wood = v.spar_wood_thickness * v.spar_wood_width * v.proj_span
		v.mass_spar_wood = v.vol_spar_wood * v.balsa_density
		v.vol_spar_cap = 2 * v.spar_cap_thickness * v.proj_span * v.spar_wood_width 
		v.mass_spar_cap = v.vol_spar_cap * v.carbon_fiber_density * v.fiberglass_layup_epoxy_factor
		v.mass_glass = (v.spar_wood_thickness + v.spar_wood_width) * v.proj_span * 2 * v.fiberglass_areal_density * v.fiberglass_layup_epoxy_factor
		v.mass_spar_misc = 0.009
		v.mass_spar = v.mass_spar_cap + v.mass_spar_wood + v.mass_glass + v.mass_spar_misc

		v.total_mass += v.mass_spar
		v.masses["Spar"] = asb.MassProperties(mass = v.mass_spar, x_cg = v.c_r / 4)

	def spar_weight_no_af(v):
		
		v.spar_cap_thickness = 0.0004 #2mm worth of carbon?
		v.spar_wood_thickness = 0.1 * v.c_r - v.spar_cap_thickness
		v.spar_wood_width = 0.0127 #1/2 inch balsa for now
		v.vol_spar_wood = v.spar_wood_thickness * v.spar_wood_width * v.proj_span
		v.mass_spar_wood = v.vol_spar_wood * v.balsa_density
		v.vol_spar_cap = 2 * v.spar_cap_thickness * v.proj_span * v.spar_wood_width 
		v.mass_spar_cap = v.vol_spar_cap * v.carbon_fiber_density * v.fiberglass_layup_epoxy_factor
		v.mass_glass = (v.spar_wood_thickness + v.spar_wood_width) * v.proj_span * 2 * v.fiberglass_areal_density * v.fiberglass_layup_epoxy_factor
		v.mass_spar_misc = 0.009
		v.mass_spar = v.mass_spar_cap + v.mass_spar_wood + v.mass_glass + v.mass_spar_misc

		v.total_mass += v.mass_spar
		v.masses["Spar"] = asb.MassProperties(mass = v.mass_spar, x_cg = v.c_r / 4)

	def tail_weight(v):
		v.vol_tail = v.hor_stab.volume() + v.vert_stab.volume()
		v.mass_foam_tail = v.vol_tail * v.foam_density

		v.area_wetted_tail = v.hor_stab.area(type="wetted") + v.vert_stab.area(type="wetted")
		v.mass_wetted_tail = v.area_wetted_tail * v.fiberglass_areal_density * v.fiberglass_layup_epoxy_factor

		v.mass_tail = v.mass_foam_tail + v.mass_wetted_tail

		v.total_mass += v.mass_tail
		v.masses["Tail"] = asb.MassProperties(mass = v.mass_tail, x_cg = v.x_tail + v.c_h / 3)

	def wing_weight_no_af(v): #not including the spar
		v.area_chord = 0.07 #from xfoil
		
		v.n_stringer = 10
		v.stringer_width = 0.007
		v.stringer_thickness = 0.001
		v.vol_stringer_wood = v.proj_span * v.n_stringer * v.stringer_width * v.stringer_thickness 

		v.n_ribs = 7 #for one half of the plane
		v.rib_chords = 2 * [v.c_r + (v.c_t - v.c_r) * i / (v.n_ribs - 1) for i in range(v.n_ribs)]
		v.vol_wood = (0.75 * v.plywood_thickness * v.area_chord * sum(c**2 for c in v.rib_chords))
		v.mass_plywood = v.plywood_density * v.vol_wood
		v.mass_stringer = v.vol_stringer_wood * v.basswood_density
		v.mass_wood = v.mass_plywood + v.mass_stringer
		v.mass_wood = v.mass_wood * 1.35 #correction factor for extra glue

		v.vol_foam = 0.2 * v.area_chord * v.b  # just the control surfaces, not that much right
		v.vol_control_wood = 0.08 * v.c_r * v.plywood_thickness * v.proj_span
		v.mass_control_wood = v.vol_control_wood * v.plywood_density
		v.mass_foam = v.vol_foam * v.foam_density + v.mass_control_wood
		v.mass_foam = v.mass_foam * 1.2

		v.skin_monokote = v.c_r * (0.7 - 0.25) * 2 * 1.1 * v.proj_span * v.monokote_density #1.1 is correction for curvature
		v.skin_carbon_fiber = v.c_r * 0.25 * 1.3 * 2 * v.proj_span * v.carbon_fiber_areal_density #1.3 is correction for curvature
		v.mass_skin = v.skin_monokote + v.skin_carbon_fiber
		v.mass_skin = v.mass_skin * 1.1 #1.1 correction factor

		v.mass_wing = v.mass_wood + v.mass_foam + v.mass_skin
		v.total_mass += v.mass_wing
		v.masses["Wing"] = asb.MassProperties(mass = v.mass_wing, x_cg = v.c_r/3) #check where c_g of wing would be

	def tail_weight_no_af(v):
		v.area_chord_tail = 0.07 #a guess for area of chord length 1

		v.vol_tail = v.area_chord_tail * v.c_h ** 2 * (v.b_h + v.b_v)
		v.mass_foam_tail = v.vol_tail * v.foam_density

		v.area_wetted_tail = v.c_v * 2.2 * v.b_v + v.c_h * 2.2 * v.b_h #this is my guess for area two sides, around chord length 1.1x per side?
		v.mass_wetted_tail = v.area_wetted_tail * v.fiberglass_areal_density * v.fiberglass_layup_epoxy_factor

		v.mass_tail = v.mass_foam_tail + v.mass_wetted_tail
		v.total_mass += v.mass_tail
		v.masses["Tail"] = asb.MassProperties(mass = v.mass_tail, x_cg = v.x_tail + v.c_h / 3)

	def weights(self):
		positions(self)
		wing_weight(self)
		misc_weight(self)
		tail_weight(self)
		boom_weight(self)
		spar_weight(self)

	def weights_no_af(self):
		positions(self)
		wing_weight_no_af(self)
		misc_weight(self)
		tail_weight_no_af(self)
		boom_weight(self)
		spar_weight_no_af(self)

#MISSION DEFINITIONS

	def M2(v):
		v.M2_score = v.mass_sensor_total / v.M2_time
		return v.M2_score	

	def M3(v): 
		v.M3_laps = v.velocity_M3 * 5 * 60 / v.lap_dist
		v.M3_score = v.mass_sensor / v.M3_laps
		return v.M3_score

	def both_missions(v):
		v.both_mission_score = M2(v) / v.M2_max + M3(v) / v.M3_max
		return v.both_mission_score

#PROPULSION VAGUE

	def run_propulsion(v):
		general_power(v)
		power_M2(v)
		power_M3(v)

	def general_power(v):
		v.battery_power = 100 * 3600 #J
		v.safety_factor = 2

	def power_M2(v):
		v.M2_distance = v.lap_dist * 5 # 5 laps
		v.M2_time = v.M2_distance / v.velocity_M2

		v.opti.subject_to(v.M2_time * v.drag_M2 * v.velocity_M2 < v.battery_power / v.safety_factor)

	def power_M3(v):
		v.opti.subject_to(v.drag_M3 * 60 * 5 * v.velocity_M3 < v.battery_power / v.safety_factor)

#SENSOR INITIALIZATION

	def sensor_dim(v):
		v.sensor_height = v.opti.variable(init_guess = 0.1, lower_bound = 3 * units.inch, upper_bound = 6 * units.inch)
		v.sensor_width = v.opti.variable(init_guess = 0.1, lower_bound = 3 * units.inch, upper_bound = 6 * units.inch)
		v.sensor_length = v.opti.variable(init_guess = 8 * units.inch, lower_bound = 6 * units.inch, upper_bound = 12 * units.inch)

		v.sensor_vol = v.sensor_length * v.sensor_width * v.sensor_height

		v.sensor_sim_number = 2

	def sensor_weight(v):
		v.ratio_container_sensor = 0.5
		v.mass_sensor = v.opti.variable(init_guess = 1, lower_bound = 0) #just the sensor
		v.mass_container = v.opti.variable(init_guess= 1, lower_bound = v.mass_sensor * v.ratio_container_sensor) #just the sensor
		v.mass_sim_sensor = (v.mass_sensor + v.mass_container) * v.sensor_sim_number
		v.mass_sensor_total = v.mass_sim_sensor + v.mass_sensor + v.mass_container

	def sensor_constraint(v):
		v.opti.subject_to((v.mass_sensor + v.mass_container) / (v.sensor_vol) < 3000) #less than straight steel

	def create_sensor(v):
		sensor_dim(v)
		sensor_weight(v)
		sensor_constraint(v)

#STRUCTURE INITIALIZATION

	def run_structures(v):
		wing_bending(v)

	def torsion(opti, params):
		#checks the torsion of the tail to the structures of the wing
		
	def wing_bending(v):
		v.spar_I = 2 * (v.spar_wood_thickness * v.spar_cap_thickness**3 / 12 + (v.spar_wood_thickness/2) **2 * v.spar_cap_thickness * v.spar_wood_width)
		v.N_max = 10
		v.max_load_per_length_M2 = v.N_max * v.g * v.mass_M2 / v.proj_span

		v.deflection_M2 = v.max_load_per_length_M2 * (v.b / 2)**4 / (8 * v.carbon_fiber_youngs_modulus * v.spar_I)
		v.deflection_max = 0.05
		v.opti.subject_to(v.deflection_M2 < v.deflection_max)

if __name__ == "__main__":
	v = variables()
	sol=v.optimize()

	for name, value in sorted(vars(v).items()):
		# print(f"{name}: {value}")
		print(f"{name}: {sol.value(value)}")