from opt import variables
import casadi as ca
import numpy 
import matplotlib.pyplot as plt


ca.GlobalOptions.setNumpyMode(-1)
S_arr= []
score_arr =[]
m2_arr = []
m3_arr = []
for S in numpy.linspace(1, 5, 30):
	
	v = variables()
	v.AR = S
	s_label = "Aspect Ratio"
	print(S)
	try:
		sol=v.optimize(verbose=False)

		# Only extract values if solve succeeded
		S_arr.append(S)
		# print(S)
		# print(sol.value(v.both_mission_score))

		score_arr.append(sol.value(v.both_mission_score))
		m2_arr.append(sol.value(v.M2_score))
		m3_arr.append(sol.value(v.M3_score))

	except RuntimeError:
		# Solve failed -> don't extract anything
		continue


import numpy as np
S_arr = np.asarray(S_arr)
score_arr = np.asarray(score_arr)
m2_arr = np.asarray(m2_arr)
m3_arr = np.asarray(m3_arr)

# Find the midpoint of the S sweep (closest point to the center value)
S_center = (S_arr.min() + S_arr.max()) / 2
mid_idx = np.argmin(np.abs(S_arr - S_center))
S_mid = S_arr[mid_idx]

print(f"Midpoint S = {S_mid} at index {mid_idx}")

def pct_diff_from_mid(arr, mid_idx):
    mid_val = arr[mid_idx]
    return (arr - mid_val) / mid_val * 100

score_pct = pct_diff_from_mid(score_arr, mid_idx)
m2_pct = pct_diff_from_mid(m2_arr, mid_idx)
m3_pct = pct_diff_from_mid(m3_arr, mid_idx)

# Quick look
# for name, pct in [("both", score_pct), ("m2", m2_pct), ("m3", m3_pct)]:
#     print(f"\n{name}: midpoint value = {pct[mid_idx]:.4f} (should be 0%)")
#     for S, p in zip(S_arr, pct):
#         print(f"  S={S:.4g}: {p:+.2f}%")

plt.plot(S_arr, score_pct, label="both")
plt.plot(S_arr, m2_pct, label="m2")
plt.plot(S_arr, m3_pct, label="m3")
plt.axhline(0, color='k', linewidth=0.5, linestyle='--')
plt.xlabel(s_label)
plt.ylabel("% difference from midpoint")
plt.title(f"% difference score vs {s_label}")
plt.legend()
plt.show()

		