from numpy import *
from matplotlib.pyplot import *
import scipy.integrate as integrate
from scipy import heaviside
import matplotlib.cm as cm
from functions_engraving_study import *

""" Plot results """
fig, ax = subplots()
dR0_values = [.005]
dT_max_values = [1000, ] # Here it is a global value, for y_s = 1 !

cms = (cm.Purples, cm.Blues, cm.Greens, cm.Oranges, cm.Reds)
colors = []
for ii in range(len(dR0_values)):
    colors.append(cms[ii](linspace(.5, 1, len(dT_max_values))))
yys = linspace(0, 1, 11)
ii = 0
for ii, dR0 in enumerate(dR0_values):
    for jj, dT_max in enumerate(dT_max_values):
        res_01 = []
        res_10 = []
        for y_s in yys:
            res_01.append(dOmega1_01(y_s, dR0, dT_max))
            res_10.append(dOmega1_10(y_s, dR0, dT_max))
        res_01 = array(res_01)
        res_10 = array(res_10)

        plot(yys, res_01, '->', c=colors[ii][jj], ms=10, label='$\Delta T_{max} = %s$\n$\delta_{\\rho_0} = %s$'%(dT_max, dR0))#'$\delta_{1^{\rightarrow}}$')
        plot(yys, res_10, '-<', c=colors[ii][jj], ms=10, label='$\Delta T_{max} = %i$\n$\delta_{R_0} = %.3f$'%(dT_max, dR0))#label='$\delta_{1^{\lefttarrow}}$')

ax.axhline(0, linestyle='--', color='k')
ax.set_xlabel('Heat point position (norm.)')
ax.set_ylabel('$\delta_{\Omega_n}$')
ax.set_title('$\Delta T_{max} = %i ;; \delta_{R_0} = %.1f$'%(dT_max, dR0))
ax.legend(ncol=len(dR0_values))

show()
