from opt import variables
import casadi as ca
import numpy 
import matplotlib.pyplot as plt


ca.GlobalOptions.setNumpyMode(-1)
S_arr= []
score_arr =[]
m2_arr = []
m3_arr = []
for S in numpy.linspace(1, 10,50):
	
	v = variables()
	v.payload_mass_frac = S
	try:
		sol=v.optimize(verbose=False)

		# Only extract values if solve succeeded
		S_arr.append(S)
		print(S)
		# print(sol.value(v.both_mission_score))

		score_arr.append(sol.value(v.both_mission_score))
		m2_arr.append(sol.value(v.M2_score))
		m3_arr.append(sol.value(v.M3_score))

	except RuntimeError:
		# Solve failed -> don't extract anything
		continue


plt.plot(S_arr, score_arr, label="both")
plt.plot(S_arr, m2_arr, label="m2")
plt.plot(S_arr, m3_arr, label="m3")
plt.legend()
plt.show()

		