import numpy as np
from scipy.optimize import minimize, Bounds
import matplotlib.pyplot as plt

def opt_obj():
	S_lower, S, S_upper = [0, 0.1, 0.2]
	AR_lower, AR, AR_upper = [0, 4, 10]
	WS_lower, WS, WS_upper = [0, 3, 5]
	TS_lower, TS, TS_upper = [0, 1, 2]

	initial_guess = (S, AR, WS, TS)
	bounds = Bounds(lb = (S_lower, AR_lower, WS_lower, TS_lower), ub = (S_upper, AR_upper, WS_upper, TS_upper))

	method = "SLSQP"

	constra

	