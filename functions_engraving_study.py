from numpy import *
from scipy import heaviside
import scipy.integrate as integrate

""" Mode shapes """
length = 1 # We need to normalize for this to work !
beta1 = pi * 0.59686 / length
C1 = (cos(beta1*length)+cosh(beta1*length))/(sin(beta1*length)+sinh(beta1*length))
beta2 = pi * 1.49418 / length
C2 = (cos(beta2*length)+cosh(beta2*length))/(sin(beta2*length)+sinh(beta2*length))
def phi1(y): # Mode shape for mode 1
    return cosh(beta1*y) - cos(beta1*y) + C1 * (sin(beta1*y)-sinh(beta1*y))
def phi1_pp(y): # Second derivative
    return beta1**2 * (cosh(beta1*y) + cos(beta1*y) - C1 * (sin(beta1*y) + sinh(beta1*y)))
def phi2(y): # Mode shape for mode 2
    return cosh(beta2*y) - cos(beta2*y) + C2 * (sin(beta2*y)-sinh(beta2*y))
def phi2_pp(y): # Second derivative
    return beta2**2 * (cosh(beta2*y) + cos(beta2*y) - C2 * (sin(beta2*y) + sinh(beta2*y)))

"""
indexes 01 : clamped end to free end (y = 0 to y = 1)
indexes 10 : free end to clamped end (y = 1 to y = 0)
"""
alphaE = -1.2e-4 # 1/E dE/dT
int1a = integrate.quad(lambda y:phi1_pp(y)**2, 0, 1)[0]
int1b = integrate.quad(lambda y:phi1(y)**2, 0, 1)[0]
int2a = integrate.quad(lambda y:phi2_pp(y)**2, 0, 1)[0]
int2b = integrate.quad(lambda y:phi2(y)**2, 0, 1)[0]


""" Direction and mode invariant : temperature """
def dT(y, y_s, dT_max):
    # Here T_max is a global value, the temperature reached when y_s = 1 !
    if y_s == 0:
        return 0
    else:
        return y/y_s * (dT_max*y_s) * heaviside(-(y - y_s), .5) + (dT_max*y_s) * heaviside(y - y_s, .5)

""" First Mechanical mode """
# Direction invariant : temperature
def dOmega1_T(y_s, dT_max):
    intT = integrate.quad(lambda y:dT(y, y_s, dT_max)*phi1_pp(y)**2, 0, 1)[0]
    return .5*alphaE*intT/int1a

# Sens base -> bout
def dR_01(y, y_s, dR0):
    return dR0 * (heaviside(y - y_s, 1) - 1)

def dOmega1_R_01(y_s, dR0):
    int_a = integrate.quad(lambda y:dR_01(y, y_s, dR0)*phi1_pp(y)**2, 0, 1)[0]
    int_b = integrate.quad(lambda y:dR_01(y, y_s, dR0)*phi1(y)**2, 0, 1)[0]
    return 2*int_a/int1a - int_b/int1b

def dOmega1_01(y_s, dR0, dT_max):
    return dOmega1_T(y_s, dT_max) + dOmega1_R_01(y_s, dR0)

# Sens bout -> base
def dR_10(y, y_s, dR0):
    return - dR0 * heaviside(y - y_s, 1)

def dOmega1_R_10(y_s, dR0):
    int_a = integrate.quad(lambda y:dR_10(y, y_s, dR0)*phi1_pp(y)**2, 0, 1)[0]
    int_b = integrate.quad(lambda y:dR_10(y, y_s, dR0)*phi1(y)**2, 0, 1)[0]
    return 2*int_a/int1a - int_b/int1b

def dOmega1_10(y_s, dR0, dT_max):
    return dOmega1_T(y_s, dT_max) + dOmega1_R_10(y_s, dR0)

""" Second mechanical mode """
# Direction invariant : temperature
def dOmega2_T(y_s, dT_max):
    intT = integrate.quad(lambda y:dT(y, y_s, dT_max)*phi2_pp(y)**2, 0, 1)[0]
    return .5*alphaE*intT/int1a

def dOmega2_R_01(y_s, dR0):
    int_a = integrate.quad(lambda y:dR_01(y, y_s, dR0)*phi2_pp(y)**2, 0, 1)[0]
    int_b = integrate.quad(lambda y:dR_01(y, y_s, dR0)*phi2(y)**2, 0, 1)[0]
    return 2*int_a/int2a - int_b/int2b
def dOmega2_01(y_s, dR0, dT_max):
    return dOmega2_T(y_s, dT_max) + dOmega2_R_01(y_s, dR0)


def dOmega2_R_10(y_s, dR0):
    int_a = integrate.quad(lambda y:dR_10(y, y_s, dR0)*phi2_pp(y)**2, 0, 1)[0]
    int_b = integrate.quad(lambda y:dR_10(y, y_s, dR0)*phi2(y)**2, 0, 1)[0]
    return 2*int_a/int2a - int_b/int2b
def dOmega2_10(y_s, dR0, dT_max):
    return dOmega2_T(y_s, dT_max) + dOmega2_R_10(y_s, dR0)


""" Density change """
# Direction invariant : temperature
def dOmega1_T(y_s, dT_max):
    intT = integrate.quad(lambda y:dT(y, y_s, dT_max)*phi1_pp(y)**2, 0, 1)[0]
    return .5*alphaE*intT/int1a

# Sens base -> bout
def drho_01(y, y_s, drho0):
    return drho0 * (heaviside(y - y_s, 1) - 1)

def dOmega1_rho_01(y_s, drho0):
    int_b = integrate.quad(lambda y:dR_01(y, y_s, drho0)*phi1(y)**2, 0, 1)[0]
    return .5 * int_b/int1b
def dOmega1_density_01(y_s, dR0, dT_max):
    return dOmega1_T(y_s, dT_max) + dOmega1_R_01(y_s, dR0)

# Sens bout -> base
def dR_10(y, y_s, dR0):
    return - dR0 * heaviside(y - y_s, 1)

def dOmega1_R_10(y_s, dR0):
    int_a = integrate.quad(lambda y:dR_10(y, y_s, dR0)*phi1_pp(y)**2, 0, 1)[0]
    int_b = integrate.quad(lambda y:dR_10(y, y_s, dR0)*phi1(y)**2, 0, 1)[0]
    return 2*int_a/int1a - int_b/int1b
def dOmega1_10(y_s, dR0, dT_max):
    return dOmega1_T(y_s, dT_max) + dOmega1_R_10(y_s, dR0)

