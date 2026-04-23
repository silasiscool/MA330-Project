from sympy import N, lambdify
from scipy.integrate import quad

# Calculate RMS error between two symbolic functions
def RMSE(f1, f2, t, P):
    f = lambdify(t, (f1-f2)**2) # Convert inner function to lambda function
    integralResult, E = quad(f, -P/2, P/2, limit=100) # Use integral approximation to find error
    return N(1/P*integralResult) # Finish RMS formula
    # return N(1/P*integrate((f1-f2)**2,(t, -P/2,P/2)))

# Calculate Normalised RMS error between two symbolic functions
def NRMSE(f1, f2, t, P, A):
    return RMSE(f1,f2,t,P) / A