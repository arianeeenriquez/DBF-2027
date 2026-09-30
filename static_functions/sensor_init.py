from aerosandbox.tools import units

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
