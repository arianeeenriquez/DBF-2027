import aerosandbox.numpy as np
from aerosandbox.tools import units

def run_aero_no_af(v):
	op_point(v)
	CL_no_af(v)
	CD_plane_no_af(v)
	drag_no_af(v
)
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
