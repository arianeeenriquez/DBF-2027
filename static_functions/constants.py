#a file of all constants

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
