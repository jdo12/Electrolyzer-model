# -*- coding: utf-8 -*-
# !/usr/bin/env python3
# PEM Electrolyzer model v0.1.4
"""
0.1.2: Includes Compressor model and energy efficiency based on input Power curves
0.1.3: Generic model that fits the user inputted power curve
0.1.4: Functional optimizations
"""
from math import log10, exp
from numpy import zeros
from scipy.optimize import fsolve
from helpers import fit_pdat, vi_calc
from matplotlib.pyplot import plot, title, show

# Model parameters to be included in config.ini
T = 25.0                                            # (°C) Electrolyzer operating temperature
Pe = 101000                                         # (pa)
Pe_out = 50000000                                   # (Pa) state of maximum charge in the tank
eta_c = 0.8                                         # compressor's isentropic efficiency
A = 0.25                                            # (m^2) area of electrode
Nc = 12                                             # Number of cells connected in series
F = 96485.34                                        # (C/mol) Faraday's constant
ne = 2 
charge_lvl_min = 20                                 # (%) minimum state of charge in the hydrogen tank
charge_lvl_max = 95                                 # (%) maximum state of charge in the hydrogen tank
soc_lvl_init = 19                                   # (%) initial state of charge
DG_25 = 237000.0
DG_80 = 228480.0
DH = 286000.0
R = 8.31445                                         # (J/mol-K) universal constant of gases
Ne = 50                                             # Number of Electrolyzers
V_tank = 0.3                                        # (m^3) volume of the tank
Nt = 3                                              # number of tanks to be charged
r1 = 7.331e-5                                       # (ohm m^2) ri parameter for ohmic resistance of electrolyte
r2 = -1.107e-7                                      # (ohm m2 °C^-1)
r3 = 0
s1 = 1.586e-1                                       # (V) si and ti parameters for over-voltage on electrodes
s2 = 1.378e-3                                       # (V°C^-1)
s3 = -1.606e-5                                      # V °C^-2)
t1 = 1.599e-2
t2 = -1.302
t3 = 4.213e2
# V vs i characteristic curve
I_initial = 0
I_final = 870
I_step = 1
a1 = 0.995                                          # 99.5 %
a2 = -9.5788                                        # (m ^ 2 * A ^ -1)
a3 = -0.0555                                        # (m ^ 2 * A ^ -1 *°C)
a4 = 0
a5 = 1502.7083                                      # (m ^ 4 * A ^ -1)
a6 = -70.8005                                       # (m ^ 4 * A ^ -1 *°C-1)
a7 = 0
gamma = 1.41
cpH2 = 14.31                                        # kj / kg - K
x0_1 = 1.6
x0_2 = 80
power_data = 'pdata.csv'                            # Power curve time-series CSV data


# Compute starting parameters for model init
min_charge = charge_lvl_min*1e-2*Pe_out             # min. state of charge in the hydrogen tank
max_charge = charge_lvl_max*1e-2*Pe_out             # max. state of charge in the hydrogen tank
soc_i = soc_lvl_init*1e-2*Pe_out                    # SOC initial
n_i = soc_i*V_tank/R/(T+273.15)
DG = DG_25-(T-25)/55*(DG_25-DG_80)
V_rev = DG/ne/F
V_init = round(V_rev, 2)
Vtn = DH/ne/F                                       # thermo-neutral voltage

nn = int((I_final-I_initial)/I_step)
I = [i for i in range(0, I_final + 1, I_step)]

# Compute V and Id
V_graph = [V_rev + (r1+r2*T)*I[i]/A+(s1+s2*T+s3*T**2)*log10((t1+t2/T+t3/T**2)*I[i]/A+1)
           for i in range(int(nn) + 1)]
Id = [i/2500 for i in I]

# Fit the power curve data
p, timespan = fit_pdat(power_data)

Pl, Pr, P_tot, Ir, P, V, ne, Qh2_V, Qh2_m,\
    m_dotH2, P_tank, Tout, W_c, soc = [zeros((timespan-1, 1), dtype=float) for i in range(14)]
moles = zeros((timespan, 1), dtype=float)
moles[0][0] = n_i

for i in range(timespan-1):
    # Power profile (kW)
    Pl[i] = sum([p[j]*(i+1)**(21-j) for j in range(22)])
    Pr[i] = 1e3*Pl[i]/Ne                            # (W) to power one stack of electrolyzer
 
    V[i], Ir[i] = fsolve(vi_calc, [x0_1, x0_2], args=(Pr[i], Nc, V_rev, T, r1, r2, s1, s2, s3, t1, t2, t3, A))
    P[i] = Ne*Nc*V[i]*Ir[i]

    # Compute Faraday Efficiency
    nf = a1*exp((a2+a3*T+a4*T**2)/(Ir[i]/A)+(a5+a6*T+a7*T**2)/(Ir[i]/A)**2)

    # Energy (or voltaje) efficiency of a cell
    ne[i] = Vtn/V[i]

    # Flow of H2 produced
    Qh2_V[i] = 80.69*Nc*Ir[i]*nf/2/F                                # (Nm^3/h)
    Qh2_m[i] = Ne*Nc*Ir[i]*nf/2/F                                   # (mol/s)
    m_dotH2[i] = Qh2_m[i]*0.001                                     # (kg/s)

    # Compressor model ERROR IN MOLES INDEX ERROR
    P_tank[i] = moles[i]*R*(T+273.15)/V_tank
    Tout[i] = (T+273.15)*(P_tank[i]/Pe)**((gamma-1)/gamma)
    W_c[i] = (m_dotH2[i]/eta_c)*cpH2*(Tout[i]-(T+273.15))           # (kW)
    P_tot[i] = W_c[i]*1000 + P[i]

    # Storage Tank
    # Number of moles in time i in the tank
    moles[i+1] = moles[i]+Qh2_m[i]*1/Nt
    soc[i] = P_tank[i]

    if soc[i] >= max_charge:
        print('Tank is fully charged!!\n'
              'Charge time:%d sec.' % i)

        soc_i = round((soc[i] / max_charge * 100)[0], 1)
        print("State of Charge:%d%%" % soc_i)
        break


title('Pl')    
plot(Pl)
show()    

title('W_c')    
plot(W_c[:i])
show()
title('Ptank')
plot(P_tank[:i])
show()
title('Ir')
plot(Ir[:i])
show()
title('V')
plot(V[:i])
show()
title('moles')
plot(moles[:i])
show()
title('m_dotH2')
plot(m_dotH2[:i]*1000)                                              # (g/s)
show()
title('Ptot')
plot(P_tot[:i])
show()
title('Pr')
plot(P[:i], 'r')
show()
