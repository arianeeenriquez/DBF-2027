import aerosandbox as asb
import aerosandbox.numpy as np

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