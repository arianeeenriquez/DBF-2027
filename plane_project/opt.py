import aerosandbox as asb
import aerosandbox.numpy as np
import json
from aerosandbox.tools import units


class variables:
	def __init__(self):
		self.opti = asb.Opti()
		self.order = 2


		#overall wing design
		self.S = self.opti.variable(init_guess = 0.2, lower_bound = 0.01, upper_bound = 1)
		self.AR = self.opti.variable(init_guess = 3, lower_bound = 2, upper_bound = 15)
		# self.AR = (5.79 * units.foot) ** 2 / self.S
		self.taper = 1
		self.airfoil = asb.Airfoil("sd7032")

		self.dihedral = 5 #degrees
		self.washout = 0 #degrees
		#tail design
		self.AR_v = self.opti.variable(init_guess = 2, lower_bound = 1, upper_bound = 5)
		self.AR_h = self.opti.variable(init_guess = 4, lower_bound = 1, upper_bound = 5)
		self.S_v = self.opti.variable(init_guess = 0.01, lower_bound = 0.001, upper_bound = .1)
		self.S_h = self.opti.variable(init_guess = 0.01, lower_bound = 0.001, upper_bound = .1)

		self.decalage = -5 #degree

		self.airfoil_tail = asb.Airfoil("naca0012")
		#locations of things
		self.x_tail =  self.opti.variable(init_guess = 1, lower_bound = 0.7, upper_bound = 1.2)

		#wing details
		self.masses = {}
		self.missions = ("M2", "M3")
		self.avl_mass = {mission: asb.MassProperties(mass=0) for mission in self.missions} #running totals for AVL mass file export, per mission
		self.total_mass = 0
		self.M2_max = 0.07
		self.M3_max = 1.14

	def _track_avl_mass(self, key, missions=None):
		#adds self.masses[key] into the running AVL mass total for each mission listed
		for mission in (missions if missions is not None else self.missions):
			self.avl_mass[mission] = self.avl_mass[mission] + self.masses[key]

	def optimize(self, verbose=True):	
		self.constants() 		#establishes known constants
		self.run_dimensions()
		self.run_aero_constraints()
		self.create_sensor()

		given_airfoil = False
		self.payload_mass_frac = self.opti.variable(init_guess = 1.5, lower_bound = 1, upper_bound = 2)

		if given_airfoil:
			self.weights()
			self.create_aero_plane()
			self.run_aero()
		else:
			self.weights_no_af()
			self.create_aero_plane()
			self.run_aero_no_af()

		self.run_propulsion()
		self.run_structures()

		self.opti.maximize(self.objective())
		sol = self.opti.solve(verbose=verbose)
		for mission in self.missions:
			self.avl_mass[mission] = sol(self.avl_mass[mission]) #convert from symbolic opti variables to solved numeric values
			self.avl_mass[mission].export_AVL_mass_file(f"example_{mission}.mass")
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

	def create_boom(self):
		self.boom = asb.Fuselage(
		    xsecs=[ 
		        asb.FuselageXSec(
		            xyz_c=[self.x_motor, 0, -units.inch/2],
		            radius=units.inch,
		        ),
		        asb.FuselageXSec(
		            xyz_c=[self.x_tail, 0,  -units.inch/2],
		            radius=units.inch,
		        ),

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
		self.create_boom()
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

	def eval_aero(self):
	    # returns an AeroBuildu

	    self.ab = {}
	    self.aero = {}
	    self.avl_analysis = {}
	    for mission in self.missions:
	        self.ab[mission] = asb.AeroBuildup(
	            airplane=self.plane,
	            op_point=self.op_pt[mission]
	        )

	        self.aero[mission] = self.ab[mission].run_with_stability_derivatives(
	            alpha=True,
	            beta=True,
	            p=False,
	            q=False,
	            r=False
	        )
	        self.avl_analysis[mission] = asb.aerodynamics.aero_3D.AVL(
	        	airplane=self.plane,
	        	op_point = self.op_pt[mission],
	        	xyz_ref = [0,0,0]
	        	)


	def lift(self):
	    #cruise: same mass build-up as CL_no_af, but lift comes straight from AeroBuildup instead of a CL formula
	    # self.mass_empty = self.total_mass
	    self.mass_empty = self.mass_sensor_total / self.payload_mass_frac

	    #payload carried differs by mission: M2 carries sensor + container, M3 carries just the sensor
	    self.payload_mass = {
	        "M2": self.mass_sensor + self.mass_container,
	        "M3": self.mass_sensor,
	    }

	    self.mass = {}
	    self.CL = {}
	    for mission in self.missions:
	        self.mass[mission] = self.mass_empty + self.payload_mass[mission]
	        self.opti.subject_to(self.mass[mission] < 60 * units.pound)
	        self.opti.subject_to(self.aero[mission]["CL"] < 0.5)
	        self.CL[mission] = self.aero[mission]["CL"]

	        #lift = weight in cruise, using AeroBuildup's own lift force instead of solving CL by hand
	        self.opti.subject_to(self.aero[mission]['L'] == self.mass[mission] * self.g)

	def drag(self):
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
		    self.drag[mission] =  self.aero[mission]['D'] + self.q[mission] * self.extra_drag_area[mission]
		    self.drag[mission] *= interference_factor

	def takeoff_constraint(self):
		self.tw_to = (1.21 / (self.g * self.rho * self.CL_max * self.S_g) * self.mass["M2"] * self.g / self.S 
			+ 0.05 ) #subject to mu and other stuff, look at 2025 design doc
		self.opti.subject_to(self.tw_to < 0.9)

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
 
	def run_aero(self):
		self.op_point()
		self.eval_aero()
		self.lift()
		self.drag()
		self.takeoff_constraint()
		self.stability_constraints()

#REDUCED ORDER AERODYNAMICS	
	def run_aero_no_af(self):
		self.op_point()
		self.CL_no_af()
		self.CD_plane_no_af()
		self.drag_no_af()
		self.takeoff_constraint()

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
		self.CL = {}
		self.N = 1 #cruise condition N = 1
		for mission in self.missions:
			self.mass[mission] = self.mass_empty + self.payload_mass[mission]
			self.opti.subject_to(self.mass[mission] < 55 * units.pound)

			self.CL[mission] = 2 * self.N * self.mass[mission] * self.g / (self.velocity[mission] ** 2 * self.S * self.rho)
			self.opti.subject_to(self.CL[mission] < 0.5) #cruise CL

	def CD_plane_no_af(self):
	    self.CDp = 0.002 
	    self.e = 0.9 #oswald efficiency

	    self.CDi = {}
	    self.CD = {}
	    for mission in self.missions:
	        self.CDi[mission] = self.CL[mission] ** 2 / (self.e * np.pi * self.AR)
	        self.CD[mission] = self.CDp + self.CDi[mission]

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
		self.x_motor = self.opti.variable(init_guess = -.3, lower_bound= -0.5, upper_bound=-0.1)
		self.x_electronics = self.opti.variable(init_guess=-0.1, lower_bound=-0.5, upper_bound=4/5*self.c_r)
		self.x_battery = self.opti.variable(init_guess=-0.1, lower_bound=-0.5, upper_bound=4/5*self.c_r) #could update to vary between M2 and M3
		self.x_battery_avionics = self.opti.variable(init_guess=-0.2, lower_bound=-0.5, upper_bound=self.c_r)

	def misc_weight(self):
		self.mass_electronics = 0.130 + .0165 + .03 + .015 #130g esc (CC 100a) + 8ch elrs rx + 30 g of wire + 10a CC BEC (2-8s)
		self.mass_battery_avionics = 0.07 #2s 1Ah lipo
		self.mass_propeller = .06 #e.g. 16-18" prop
		self.mass_tail_servo = 2*.035  #2x wing servos for tail
		self.mass_wing_servo = 4 * .035 #4x main wing servos
		self.mass_landing_gear = 0.2
		self.mass_misc = 0.4

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
		self.area_chord = 6.43354654E-02 #from xfoil
		
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

		self.vol_foam = 0.2 * self.wing.volume()  # just the control surfaces, not that much right
		self.vol_control_wood = self.airfoil.local_thickness(x_over_c = 0.6) * self.c_r * self.plywood_thickness * self.proj_span
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

	def boom_weight(self):
		self.boom_weight = (self.x_tail - self.x_motor + self.c_h / 4) / 1.27 * 0.173 #scaling from last year
		self.total_mass += self.boom_weight
		self.masses["Boom"] = asb.MassProperties(mass = self.boom_weight, x_cg = (self.x_tail + self.c_h/4 + self.x_motor) / 2 )
		self._track_avl_mass("Boom")

	def spar_weight(self):
		
		self.spar_cap_thickness = 0.0004 #2mm worth of carbon?
		self.spar_wood_thickness = self.airfoil.local_thickness(x_over_c=0.25) * self.c_r - self.spar_cap_thickness
		self.spar_wood_width = 0.0127 #1/2 inch balsa for now
		self.vol_spar_wood = self.spar_wood_thickness * self.spar_wood_width * self.proj_span
		self.mass_spar_wood = self.vol_spar_wood * self.balsa_density
		self.vol_spar_cap = 2 * self.spar_cap_thickness * self.proj_span * self.spar_wood_width 
		self.mass_spar_cap = self.vol_spar_cap * self.carbon_fiber_density * self.fiberglass_layup_epoxy_factor
		self.mass_glass = (self.spar_wood_thickness + self.spar_wood_width) * self.proj_span * 2 * self.fiberglass_areal_density * self.fiberglass_layup_epoxy_factor
		self.mass_spar_misc = 0.009
		self.mass_spar = self.mass_spar_cap + self.mass_spar_wood + self.mass_glass + self.mass_spar_misc

		self.total_mass += self.mass_spar
		self.masses["Spar"] = asb.MassProperties(mass = self.mass_spar, x_cg = self.c_r / 4)
		self._track_avl_mass("Spar")

	def spar_weight_no_af(self):
		
		self.spar_cap_thickness = 0.002 #2mm worth of carbon?
		self.spar_wood_thickness = 0.1 * self.c_r - self.spar_cap_thickness
		self.spar_wood_width = 0.0127 #1/2 inch balsa for now
		self.vol_spar_wood = self.spar_wood_thickness * self.spar_wood_width * self.proj_span
		self.mass_spar_wood = self.vol_spar_wood * self.balsa_density
		self.vol_spar_cap = 2 * self.spar_cap_thickness * self.proj_span * self.spar_wood_width 
		self.mass_spar_cap = self.vol_spar_cap * self.carbon_fiber_density * self.fiberglass_layup_epoxy_factor
		self.mass_glass = (self.spar_wood_thickness + self.spar_wood_width) * self.proj_span * 2 * self.fiberglass_areal_density * self.fiberglass_layup_epoxy_factor
		self.mass_spar_misc = 0.009
		self.mass_spar = self.mass_spar_cap + self.mass_spar_wood + self.mass_glass + self.mass_spar_misc

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
		self.area_chord = 0.07 #from xfoil
		
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

		self.vol_foam = 0.2 * self.area_chord * self.b  # just the control surfaces, not that much right
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
		self.positions()
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
		return self.M2_score	

	def M3(self): 
		self.M3_laps = self.velocity["M3"] * 5 * 60 / self.lap_dist
		self.M3_score = self.mass_sensor / self.M3_laps
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
		self.duration = {
			"M2": self.M2_distance / self.velocity["M2"],
			"M3": 5 * 60,
		}

		for mission in self.missions:
			self.opti.subject_to(self.drag[mission] * self.velocity[mission] * self.duration[mission] < self.battery_power / self.safety_factor)
		
	def general_power(self):
		self.battery_power = 100 * 3600 #J
		self.safety_factor = 2

#SENSOR INITIALIZATION

	def sensor_dim(self):
		self.sensor_height = self.opti.variable(init_guess = 0.1, lower_bound = 3 * units.inch, upper_bound = 6 * units.inch)
		self.sensor_width = self.opti.variable(init_guess = 0.1, lower_bound = 3 * units.inch, upper_bound = 6 * units.inch)
		self.sensor_length = self.opti.variable(init_guess = 8 * units.inch, lower_bound = 6 * units.inch, upper_bound = 12 * units.inch)

		self.sensor_vol = self.sensor_length * self.sensor_width * self.sensor_height


		self.fuselage = asb.Fuselage(
		    xsecs=[
		        asb.FuselageXSec(
		            xyz_c=[-0.06, 0,  -self.sensor_height/1.2],
		            radius=0,
		        ),
		        asb.FuselageXSec(
		            xyz_c=[-0.04, 0,  -self.sensor_height/1.2],
		            radius=self.sensor_height/2,
		        ),
		        asb.FuselageXSec(
		            xyz_c=[0, 0, -self.sensor_height/1.2],
		            radius=self.sensor_height/1.2,
		        ),
		        asb.FuselageXSec(
		            xyz_c=[self.sensor_length, 0, -self.sensor_height/1.2],
		            radius=self.sensor_height/1.2,
		        ),
		        asb.FuselageXSec(
		            xyz_c=[self.sensor_length + 0.04, 0,  -self.sensor_height/1.2],
		            radius=self.sensor_height/2,
		        ),
		        asb.FuselageXSec(
		            xyz_c=[self.sensor_length + 0.06, 0,  -self.sensor_height/1.2],
		            radius=0,
		        ),

		    ]
		)

	def sensor_weight(self):
		self.ratio_container_sensor = 0.5
		self.mass_sensor = self.opti.variable(init_guess = 1.5, lower_bound = 1, upper_bound = 2) #just the sensor
		self.mass_container = self.opti.variable(init_guess= 1, lower_bound = self.mass_sensor * self.ratio_container_sensor) #just the sensor
		self.mass_sensor_total = self.mass_sensor + self.mass_container

		self.x_sensor = self.opti.variable(init_guess=0, lower_bound = -0.1) #TODO: set this to the actual sensor/container CG location, this is just a placeholder
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
		self.N_max = 10
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

if __name__ == "__main__":
	import argparse
	import casadi as ca

	parser = argparse.ArgumentParser(description="A simple CLI tool example.")
	parser.add_argument("-m", "--mission", type=int, default=0, help="mission number")
	parser.add_argument("-p", "--print", type=bool, default=False, help="print")
	args = parser.parse_args()
	mission = args.mission



	ca.GlobalOptions.setNumpyMode(-1)
	v = variables()

	v.M2_max = 0.1
	v.M3_max = 0.56

	v.mission_no = mission
	sol=v.optimize(verbose=False)
	# v.print_summary(sol)
	# v.save_solution(sol, "solution.json")

	print(sol.value(v.M2_score))
	print(sol.value(v.M3_score))
	print(sol.value(v.both_mission_score))

	    # #print everything
	if args.print:
		for name, value in sorted(vars(v).items()):
			# print(f"{name}: {value}")
			print(f"{name}: {sol.value(value)}")

		sol.value(v.plane).draw_three_view()

		# sol.value(v.avl_analysis["M2"]).write_avl(filepath="M2.avl")