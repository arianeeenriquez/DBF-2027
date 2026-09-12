def run_structures(v):
	wing_bending(v)

def torsion(opti, params):
	#checks the torsion of the tail to the structures of the wing
	pass

def wing_bending(v):
	v.spar_I = 2 * (v.spar_wood_thickness * v.spar_cap_thickness**3 / 12 + (v.spar_wood_thickness/2) **2 * v.spar_cap_thickness * v.spar_wood_width)
	v.N_max = 10
	v.max_load_per_length_M2 = v.N_max * v.g * v.mass_M2 / v.proj_span

	v.deflection_M2 = v.max_load_per_length_M2 * (v.b / 2)**4 / (8 * v.carbon_fiber_youngs_modulus * v.spar_I)
	v.deflection_max = 0.05
	v.opti.subject_to(v.deflection_M2 < v.deflection_max)
	pass

def boom_bending(opti, params):
	#checks the boom deflection maybe necessary? how much does tail make the boom bend maybe
	pass

def landing_gear(opti, params):
	#maybe checks how much force is expected of the landing gear
	pass