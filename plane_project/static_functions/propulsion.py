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
