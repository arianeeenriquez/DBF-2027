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