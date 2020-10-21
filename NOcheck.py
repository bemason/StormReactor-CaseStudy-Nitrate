# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   Brooke Mason
# @Last Modified time: 2020-10-21 09:36:12

# Import required modules
from pyswmm import Simulation, Nodes, Links
from wq_toolbox.nodes import Node_Quality
from scipy.integrate import ode 
import numpy as np
import matplotlib.pyplot as plt

# Nitrate 3 CSTRs in Series
def CSTR_tank(t, C, Qin, Cin, Qout, V, k):
    dCdt = (Qin*Cin - Qout*C)/V - k*C
    return dCdt

# DO Variables
k_DO = 0.000075 # rate/5 sec
Co_DO = 10.0    # mg/L
k_ni = 0.000029*3
#----------------------------------------------------------------------#

# Uncontrolled Simulation
# Lists to store results
Wetland_inflow = []
Wetland_conc = []
Wetland_depth = []
Wetland_volume = []
Wetland_outflow = []
Wetland_cumload = []
Wetland_Cin = []
Wetland_k = []
Wetland_DO1 = []
Wetland_DO2 = []
Wetland_DO3 = []   

# Setup toolbox simulation
with Simulation("./NO_check.inp") as sim:
    # Get asset information
    Wetland = Nodes(sim)["93-49759"]
    
    # Setup dt calculation        
    start_time = sim.start_time
    last_timestep = start_time

    # Setup CSTR solver
    solver1 = ode(CSTR_tank)
    solver1.set_integrator("dopri5")
    solver2 = ode(CSTR_tank)
    solver2.set_integrator("dopri5")
    solver3 = ode(CSTR_tank)
    solver3.set_integrator("dopri5")

    # Tracking time for DO reaction
    t1 = 0
    t2 = 0
    t3 = 0

    # Step through the simulation    
    for index,step in enumerate(sim):

        # Calculate dt
        current_step = sim.current_time
        dt = (current_step - last_timestep).total_seconds()
        last_timestep = current_step

        # Get NO conc for each asset        
        Wt_p = Wetland.pollut_quality['NO']
        Wetland_conc.append(Wt_p)

        Wt_if = Wetland.total_inflow
        Wetland_inflow.append(Wt_if)
        Wt_d = Wetland.depth
        Wetland_depth.append(Wt_d)
        Wt_v = Wetland.volume
        Wetland_volume.append(Wt_v)
        Wt_of = Wetland.total_outflow
        Wetland_outflow.append(Wt_of)

        
        # Calculate DO concentration in tank layers
        # reset DO if tank is empty
        if Wt_d <= 0.01:
            t1 = 0
            t2 = 0
            t3 = 0
            DO1 = 0.0
            DO2 = 0.0
            DO3 = 0.0
            Wetland_DO1.append(DO1)
            Wetland_DO2.append(DO2)
            Wetland_DO3.append(DO3)
            k_ni = 0.0
        
        elif 0.01 < Wt_d <= 3.00:
            # Calculate DO concentration in first layer
            t1 += dt
            t2 = 0
            t3 = 0
            DO1 = Co_DO*np.exp(-k_DO*t1)
            DO2 = 0.0
            DO3 = 0.0
            # Calculate nitate reaction rate based on DO concentration
            if DO1 <= 1.0:
                k_ni = 0.000029  # [1/5 sec]
            else:
                k_ni = 0.0
            Wetland_DO1.append(DO1)
            Wetland_DO2.append(DO2)
            Wetland_DO3.append(DO3)
        
        elif 3.0 < Wt_d <= 6.0:
            # Calculate DO concentration in first two layers
            t1 += dt
            t2 += dt
            t3 = 0
            DO1 = Co_DO*np.exp(-k_DO*t1)
            DO2 = Co_DO*np.exp(-k_DO*t2)
            DO3 = 0.0
            # Calculate nitate reaction rate based on DO concentration
            if DO1 <= 1.0:
                k_ni1 = 0.000029
            else:
                k_ni1 = 0.0
            if DO2 <= 1.0:
                k_ni2 = 0.000029
            else:
                k_ni2 = 0.0
            Wetland_DO1.append(DO1)
            Wetland_DO2.append(DO2)
            Wetland_DO3.append(DO3)
        else:
            # Calculate DO concentration in all three layers
            t1 += dt
            t2 += dt
            t3 += dt
            DO1 = Co_DO*np.exp(-k_DO*t1)
            DO2 = Co_DO*np.exp(-k_DO*t2)
            DO3 = Co_DO*np.exp(-k_DO*t3)
            # Calculate nitate reaction rate based on DO concentration
            if DO1 <= 1.0:
                k_ni1 = 0.000029
            else:
                k_ni1 = 0.0
            if DO2 <= 1.0:
                k_ni2 = 0.000029
            else:
                k_ni2 = 0.0
            if DO3 <= 1.0:
                k_ni3 = 0.000029
            else:
                k_ni3 = 0.0
            k_ni = k_ni1 + k_ni2 + k_ni3
            Wetland_DO1.append(DO1)
            Wetland_DO2.append(DO2)
            Wetland_DO3.append(DO3)
        Wetland_k.append(k_ni)
        

        # Calculate NO concentration in tanks
        # Get parameters to calculate NO
        Cin=sim._model.getNodeCin("93-49759",0)
        Wetland_Cin.append(Cin)

        #Solve ODE
        if index == 0:
            solver1.set_initial_value(0.0, 0.0)
            solver1.set_f_params(Wt_if,Cin,Wt_of,Wt_v,k_ni)
            solver1.integrate(solver1.t+dt)
            solver2.set_initial_value(0.0, 0.0)
            solver2.set_f_params(Wt_if,Cin,Wt_of,Wt_v,k_ni)
            solver2.integrate(solver2.t+dt)
            solver3.set_initial_value(0.0, 0.0)
            solver3.set_f_params(Wt_if,Cin,Wt_of,Wt_v,k_ni)
            solver3.integrate(solver3.t+dt)
        else:
            solver1.set_initial_value(solver1.y, solver1.t)
            solver1.set_f_params(Wt_if,Cin,Wt_of,Wt_v,k_ni)
            solver1.integrate(solver1.t+dt)
            solver2.set_initial_value(solver2.y, solver2.t)
            solver2.set_f_params(Wt_if,solver1.y,Wt_of,Wt_v,k_ni)
            solver2.integrate(solver2.t+dt)
            solver3.set_initial_value(solver3.y, solver3.t)
            solver3.set_f_params(Wt_if,solver2.y,Wt_of,Wt_v,k_ni)
            solver3.integrate(solver3.t+dt)
        
        # Set new concentration
        sim._model.setNodePollutant("93-49759", 0, solver3.y[0])

    sim._model.swmm_end()
    print(sim.runoff_error)
    print(sim.flow_routing_error)
    print(sim.quality_error)

# Convert inflow rate from cfs to m3/s
conv_cfs_cms = [0.02832]*len(Wetland_inflow)
Wetland_inflow_m = [a*b for a,b in zip(Wetland_inflow,conv_cfs_cms)]

# Convert outflow rate from cfs to m3/s
Wetland_outflow_m = [a*b for a,b in zip(Wetland_outflow,conv_cfs_cms)]

# Calculate load each timestep
conv_mgs_kgs = [0.000001]*len(Wetland_inflow)
timestep = [5]*len(Wetland_inflow)
Wetland_load = [a*b*c*d*e for a,b,c,d,e in zip(Wetland_conc,Wetland_outflow,conv_cfs_cms, conv_mgs_kgs,timestep)]

# Calculate cumulative load
Wetland_cumload = np.cumsum(Wetland_load)

#----------------------------------------------------------------------#
# Confirm CSTR matches steady state equilibrium concentration
# Found average V,Q,r,Co for simulation above 
V = np.mean(Wetland_volume[34560:51840])
Q = np.mean(Wetland_inflow[34560:51840])
k = np.mean(Wetland_k[34560:51840])
Co = np.mean(Wetland_Cin[34560:51840])
Cn = Co /((1+(k*V/Q))**3)
print("SS Conc:", Cn)
print("Model Conc:", Wetland_conc[-1])

#----------------------------------------------------------------------#
# Data for flooding line
Wetland_concSS  = [Cn]*len(Wetland_conc)

# Plot Result
plt.plot(Wetland_conc, color='#B08CA1', linewidth=2, label='Model')
plt.plot(Wetland_concSS, 'k--', linewidth=2, label='SS Conc')
plt.ylabel("Concentration (mg/L)")
plt.xlabel("Time (days)")
plt.xticks([0,17280,34560],["0","1","2"])
plt.legend()
plt.show()


