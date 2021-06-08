# -*- coding: utf-8 -*-
# @Author: Brooke Mason
# @Date:   2020-01-15 09:57:05
# @Last Modified by:   HP
# @Last Modified time: 2021-06-08 10:54:14

# Import required modules
from pyswmm import Simulation, Nodes, Links
from scipy.integrate import ode 
import numpy as np
import matplotlib.pyplot as plt

# 3 CSTRs in Series
def CSTR_tank(t, C, Qin, Cin, Qout, V, k):
    dCdt = (Qin*Cin - Qout*C)/V - k*C
    return dCdt

# DO Variables
k_DO = 0.000278  # rate/5 sec
Cin_DO = 9.6     # mg/L

#----------------------------------------------------------------------#

# Uncontrolled Simulation
# Lists to store results
RBasin_inflow = []
RBasin_conc = []
RBasin_depth = []
#RBasin_flooding = []
RBasin_outflow = []
RBasin_cumload = []

DBasin_inflow = []
DBasin_conc = []
DBasin_depth = []
DBasin_outflow = []
DBasin_cumload = []

Wetland_inflow = []
Wetland_conc = []
Wetland_depth = []
#Wetland_flooding= []
Wetland_volume = []
Wetland_outflow = []
Wetland_cumload = []
Wetland_DO = []  
#Wtlnd_bp_inflows = []

Channel_flow = []
Channel_conc = []
Channel_depth = []
#Channel_flooding = []
Channel_cumload = []

# Setup toolbox simulation
with Simulation("./NO.inp") as sim:
    # Get asset information
    RBasin = Nodes(sim)["93-50408"]
    DBasin = Nodes(sim)["93-50404"]
    Wetland = Nodes(sim)["93-49759"]
    Wtlnd_bypass = Links(sim)["95-70294"]
    Channel = Links(sim)["95-70277"]
    
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
    solver4 = ode(CSTR_tank)
    solver4.set_integrator("dopri5")
    solver5 = ode(CSTR_tank)
    solver5.set_integrator("dopri5")
    solver6 = ode(CSTR_tank)
    solver6.set_integrator("dopri5")


    # Step through the simulation    
    for index,step in enumerate(sim):

        # Calculate dt
        current_step = sim.current_time
        dt = (current_step - last_timestep).total_seconds()
        last_timestep = current_step

        # Get NO conc for each asset        
        RB_p = RBasin.pollut_quality['NO']
        RBasin_conc.append(RB_p)
        DB_p = DBasin.pollut_quality['NO']
        DBasin_conc.append(DB_p)
        Wt_p = Wetland.pollut_quality['NO']
        Wetland_conc.append(Wt_p)
        Ch_p = Channel.pollut_quality['NO']
        Channel_conc.append(Ch_p)

        # Get flow for  each asset
        RB_if = RBasin.total_inflow
        RBasin_inflow.append(RB_if)
        RB_d = RBasin.depth
        RBasin_depth.append(RB_d)
        #RB_fl = RBasin.flooding
        #RBasin_flooding.append(RB_fl)
        RB_of = RBasin.total_outflow
        RBasin_outflow.append(RB_of)
        DB_if = DBasin.total_inflow
        DBasin_inflow.append(DB_if)
        DB_d = DBasin.depth
        DBasin_depth.append(DB_d)
        DB_of = DBasin.total_outflow
        DBasin_outflow.append(DB_of)
        Wt_if = Wetland.total_inflow
        Wetland_inflow.append(Wt_if)
        Wt_d = Wetland.depth
        Wetland_depth.append(Wt_d)
        #Wt_fl = Wetland.flooding
        #Wetland_flooding.append(Wt_fl)
        Wt_v = Wetland.volume
        Wetland_volume.append(Wt_v)
        Wt_of = Wetland.total_outflow
        Wetland_outflow.append(Wt_of)
        #Wt_bp = Wtlnd_bypass.flow
        #Wtlnd_bp_inflows.append(Wt_bp)
        Ch_f = Channel.flow
        Channel_flow.append(Ch_f)
        Ch_d = Channel.depth
        Channel_depth.append(Ch_d)

        # Calculate DO concentration in tank
        if index == 0:
            solver4.set_initial_value(0.0, 0.0)
            solver4.set_f_params(Wt_if,Cin_DO,Wt_of,Wt_v,k_DO)
            solver4.integrate(solver4.t+dt)
            solver5.set_initial_value(0.0, 0.0)
            solver5.set_f_params(Wt_if,Cin_DO,Wt_of,Wt_v,k_DO)
            solver5.integrate(solver5.t+dt)
            solver6.set_initial_value(0.0, 0.0)
            solver6.set_f_params(Wt_if,Cin_DO,Wt_of,Wt_v,k_DO)
            solver6.integrate(solver6.t+dt)
        else:
            solver4.set_initial_value(solver4.y, solver4.t)
            solver4.set_f_params(Wt_if,Cin_DO,Wt_of,Wt_v,k_DO)
            solver4.integrate(solver4.t+dt)
            solver5.set_initial_value(solver5.y, solver5.t)
            solver5.set_f_params(Wt_if,solver4.y,Wt_of,Wt_v,k_DO)
            solver5.integrate(solver5.t+dt)
            solver6.set_initial_value(solver6.y, solver6.t)
            solver6.set_f_params(Wt_if,solver5.y,Wt_of,Wt_v,k_DO)
            solver6.integrate(solver6.t+dt)

        # Save DO concentration
        DO = solver6.y[0]
        Wetland_DO.append(DO)

        if DO <= 1.0:
            k_ni = 0.000087
        else:
            k_ni = 0.0
        
        # Calculate NO concentration in tanks
        # Get parameters to calculate NO
        Cin_NO = sim._model.getNodeCin("93-49759",0)

        #Solve Nitrate ODE
        if index == 0:
            solver1.set_initial_value(0.0, 0.0)
            solver1.set_f_params(Wt_if,Cin_NO,Wt_of,Wt_v,k_ni)
            solver1.integrate(solver1.t+dt)
            solver2.set_initial_value(0.0, 0.0)
            solver2.set_f_params(Wt_if,Cin_NO,Wt_of,Wt_v,k_ni)
            solver2.integrate(solver2.t+dt)
            solver3.set_initial_value(0.0, 0.0)
            solver3.set_f_params(Wt_if,Cin_NO,Wt_of,Wt_v,k_ni)
            solver3.integrate(solver3.t+dt)
        else:
            solver1.set_initial_value(solver1.y, solver1.t)
            solver1.set_f_params(Wt_if,Cin_NO,Wt_of,Wt_v,k_ni)
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
conv_cfs_cms = [0.02832]*len(RBasin_inflow)
RBasin_inflow_m = [a*b for a,b in zip(RBasin_inflow,conv_cfs_cms)]
DBasin_inflow_m = [a*b for a,b in zip(DBasin_inflow,conv_cfs_cms)]
Wetland_inflow_m = [a*b for a,b in zip(Wetland_inflow,conv_cfs_cms)]
Channel_flow_m = [a*b for a,b in zip(Channel_flow,conv_cfs_cms)]

# Convert outflow rate from cfs to m3/s
conv_cfs_cms = [0.02832]*len(RBasin_inflow)
RBasin_outflow_m = [a*b for a,b in zip(RBasin_outflow,conv_cfs_cms)]
DBasin_outflow_m = [a*b for a,b in zip(DBasin_outflow,conv_cfs_cms)]
Wetland_outflow_m = [a*b for a,b in zip(Wetland_outflow,conv_cfs_cms)]

# Convert flooding rate from cfs to m3/s
#conv_cfs_cms = [0.02832]*len(RBasin_inflow)
#RBasin_flooding_m = [a*b for a,b in zip(RBasin_flooding,conv_cfs_cms)]
#Wetland_flooding_m = [a*b for a,b in zip(Wetland_flooding,conv_cfs_cms)]
#Wtlnd_bypass_m = [a*b for a,b in zip(Wtlnd_bp_inflows,conv_cfs_cms)]

# Convert depth from ft to m
conv_ft_m = [0.3048]*len(RBasin_inflow)
RBasin_depth_m = [a*b for a,b in zip(RBasin_depth,conv_ft_m)]
DBasin_depth_m = [a*b for a,b in zip(DBasin_depth,conv_ft_m)]
Wetland_depth_m = [a*b for a,b in zip(Wetland_depth,conv_ft_m)]
Channel_depth_m = [a*b for a,b in zip(Channel_depth,conv_ft_m)]

# Calculate load each timestep
conv_mgs_kgs = [0.000001]*len(RBasin_inflow)
timestep = [5]*len(RBasin_inflow)
RBasin_load = [a*b*c*d*e for a,b,c,d,e in zip(RBasin_conc,RBasin_outflow,conv_cfs_cms, conv_mgs_kgs,timestep)]
DBasin_load = [a*b*c*d*e for a,b,c,d,e in zip(DBasin_conc,DBasin_outflow,conv_cfs_cms,conv_mgs_kgs,timestep)]
Wetland_load = [a*b*c*d*e for a,b,c,d,e in zip(Wetland_conc,Wetland_outflow,conv_cfs_cms, conv_mgs_kgs,timestep)]
Channel_load = [a*b*c*d*e for a,b,c,d,e in zip(Channel_conc,Channel_flow,conv_cfs_cms,conv_mgs_kgs,timestep)]

# Calculate cumulative load (dt = 1)
RBasin_cumload = np.cumsum(RBasin_load)
DBasin_cumload = np.cumsum(DBasin_load)
Wetland_cumload = np.cumsum(Wetland_load)
Channel_cumload = np.cumsum(Channel_load)

#----------------------------------------------------------------------#
# Controlled Simulation 
# Lists to store results
RBasin_inflowC = []
RBasin_concC = []
RBasin_depthC = []
#RBasin_floodingC = []
RBasin_valveC = []
RBasin_outflowC = []
RBasin_cumloadC = []

DBasin_inflowC = []
DBasin_concC = []
DBasin_depthC = []
DBasin_outflowC = []
DBasin_cumloadC = []

Wetland_inflowC = []
Wetland_concC = []
Wetland_depthC = []
#Wetland_floodingC = []
Wetland_volumeC = []
Wetland_valveC = []
Wetland_outflowC = []
Wetland_cumloadC = []
Wetland_DOC = []
#Wtlnd_bp_inflowsC = []

Channel_flowC = []
Channel_concC = []
Channel_depthC = []
Channel_outflowC = []
Channel_cumloadC = []

# DO Variables
k_DO = 0.000278  # rate/5 sec
Cin_DO = 9.6     # mg/L

# Setup toolbox simulation
with Simulation("./NO.inp") as sim:
    
    # Get asset information
    RBasin = Nodes(sim)["93-50408"]
    RB_valve = Links(sim)["95-70951"]
    DBasin = Nodes(sim)["93-50404"]
    Channel = Links(sim)["95-70277"]
    Wetland = Nodes(sim)["93-49759"]
    Wtlnd_bypass = Links(sim)["95-70294"]
    Wtlnd_valve = Links(sim)["95-70293"]

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
    solver4 = ode(CSTR_tank)
    solver4.set_integrator("dopri5")
    solver5 = ode(CSTR_tank)
    solver5.set_integrator("dopri5")
    solver6 = ode(CSTR_tank)
    solver6.set_integrator("dopri5")

    # Tracking time for control actions every 15 minutes (5 sec time step)
    _tempcount = 180

    # Step through the simulation    
    for index,step in enumerate(sim):

        # Calculate dt
        current_step = sim.current_time
        dt = (current_step - last_timestep).total_seconds()
        last_timestep = current_step

        # Get NO conc for each asset        
        RB_p = RBasin.pollut_quality['NO']
        RBasin_concC.append(RB_p)
        DB_p = DBasin.pollut_quality['NO']
        DBasin_concC.append(DB_p)
        Wt_p = Wetland.pollut_quality['NO']
        Wetland_concC.append(Wt_p)
        Ch_p = Channel.pollut_quality['NO']
        Channel_concC.append(Ch_p)

        # Get flow for  each asset
        RB_if = RBasin.total_inflow
        RBasin_inflowC.append(RB_if)
        RB_d = RBasin.depth
        RBasin_depthC.append(RB_d)
        #RB_fl = RBasin.flooding
        #RBasin_floodingC.append(RB_fl)
        RB_of = RBasin.total_outflow
        RBasin_outflowC.append(RB_of)
        DB_if = DBasin.total_inflow
        DBasin_inflowC.append(DB_if)
        DB_d = DBasin.depth
        DBasin_depthC.append(DB_d)
        DB_of = DBasin.total_outflow
        DBasin_outflowC.append(DB_of)
        Wt_if = Wetland.total_inflow
        Wetland_inflowC.append(Wt_if)
        Wt_d = Wetland.depth
        Wetland_depthC.append(Wt_d)
        #Wt_fl = Wetland.flooding
        #Wetland_floodingC.append(Wt_fl)
        Wt_v = Wetland.volume
        Wetland_volumeC.append(Wt_v)
        Wt_of = Wetland.total_outflow
        Wetland_outflowC.append(Wt_of)
        #Wt_bp = Wtlnd_bypass.flow
        #Wtlnd_bp_inflowsC.append(Wt_bp)
        Ch_f = Channel.flow
        Channel_flowC.append(Ch_f)
        Ch_d = Channel.depth
        Channel_depthC.append(Ch_d)

        # Calculate DO concentration in tank
        if index == 0:
            solver4.set_initial_value(0.0, 0.0)
            solver4.set_f_params(Wt_if,Cin_DO,Wt_of,Wt_v,k_DO)
            solver4.integrate(solver4.t+dt)
            solver5.set_initial_value(0.0, 0.0)
            solver5.set_f_params(Wt_if,Cin_DO,Wt_of,Wt_v,k_DO)
            solver5.integrate(solver5.t+dt)
            solver6.set_initial_value(0.0, 0.0)
            solver6.set_f_params(Wt_if,Cin_DO,Wt_of,Wt_v,k_DO)
            solver6.integrate(solver6.t+dt)
        else:
            solver4.set_initial_value(solver4.y, solver4.t)
            solver4.set_f_params(Wt_if,Cin_DO,Wt_of,Wt_v,k_DO)
            solver4.integrate(solver4.t+dt)
            solver5.set_initial_value(solver5.y, solver5.t)
            solver5.set_f_params(Wt_if,solver4.y,Wt_of,Wt_v,k_DO)
            solver5.integrate(solver5.t+dt)
            solver6.set_initial_value(solver6.y, solver6.t)
            solver6.set_f_params(Wt_if,solver5.y,Wt_of,Wt_v,k_DO)
            solver6.integrate(solver6.t+dt)

        # Save DO concentration
        DO_C = solver6.y[0]
        Wetland_DOC.append(DO_C)

        if DO_C <= 1.0:
            k_ni = 0.000087
        else:
            k_ni = 0.0
        
        # Calculate NO concentration in tanks
        # Get parameters to calculate NO
        Cin_NO=sim._model.getNodeCin("93-49759",0)

        #Solve ODE
        if index == 0:
            solver1.set_initial_value(0.0, 0.0)
            solver1.set_f_params(Wt_if,Cin_NO,Wt_of,Wt_v,k_ni)
            solver1.integrate(solver1.t+dt)
            solver2.set_initial_value(0.0, 0.0)
            solver2.set_f_params(Wt_if,Cin_NO,Wt_of,Wt_v,k_ni)
            solver2.integrate(solver2.t+dt)
            solver3.set_initial_value(0.0, 0.0)
            solver3.set_f_params(Wt_if,Cin_NO,Wt_of,Wt_v,k_ni)
            solver3.integrate(solver3.t+dt)
        else:
            solver1.set_initial_value(solver1.y, solver1.t)
            solver1.set_f_params(Wt_if,Cin_NO,Wt_of,Wt_v,k_ni)
            solver1.integrate(solver1.t+dt)
            solver2.set_initial_value(solver2.y, solver2.t)
            solver2.set_f_params(Wt_if,solver1.y,Wt_of,Wt_v,k_ni)
            solver2.integrate(solver2.t+dt)
            solver3.set_initial_value(solver3.y, solver3.t)
            solver3.set_f_params(Wt_if,solver2.y,Wt_of,Wt_v,k_ni)
            solver3.integrate(solver3.t+dt)
        
        # Set new concentration
        sim._model.setNodePollutant("93-49759", 0, solver3.y[0])

        # Wetland & Retention basin Control Actions (every 15 mins - 5 sec timesteps)
        if _tempcount == 180:
            # If DO level is not anoxic
            if DO_C > 1.0:
                # And if the wetland has capacity
                if Wt_d <= 9.5:
                    # Close the wetland valve and proportionally open retention basin valve C = Qmax/(A*sqrt(2*g*d))
                    Wtlnd_valve.target_setting = 0.0
                    RB_valve.target_setting = 1.75*(70.6/(np.sqrt(2*32.2*RB_d)*78.5))
                else:
                    # If not, open the wetland valve and close the RBasin valve
                    Wtlnd_valve.target_setting = 1.75*(70.6/(np.sqrt(2*32.2*Wt_d)*12.6))
                    RB_valve.target_setting = 0.0
            # If DO level is anoxic
            elif DO_C <= 1.0:
                # And if the wetland NO concentration is low, open both valves proportionally
                if solver3.y[0] <= 5.0:
                    Wtlnd_valve.target_setting = 1.75*(70.6/(np.sqrt(2*32.2*Wt_d)*12.6))
                    RB_valve.target_setting = 1.75*(70.6/(np.sqrt(2*32.2*RB_d)*78.5))
                # Else if the wetland NO concentration is high
                else:
                    # And if the wetland still has capacity, close both valves
                    if Wt_d <= 9.5:
                        Wtlnd_valve.target_setting = 0.0
                        RB_valve.target_setting = 0.0
                    # If not, open the wetland valve propotionally and close retention basin
                    else:
                        Wtlnd_valve.target_setting = 1.75*(70.6/(np.sqrt(2*32.2*Wt_d)*12.6))
                        RB_valve.target_setting = 0.0

            _tempcount= 0
        _tempcount+= 1

        Wetland_valveC.append(Wtlnd_valve.target_setting)
        RBasin_valveC.append(RB_valve.target_setting)

    sim._model.swmm_end()
    print(sim.runoff_error)
    print(sim.flow_routing_error)
    print(sim.quality_error)


# Convert inflow rate from cfs to m3/s
conv_cfs_cms = [0.02832]*len(RBasin_inflowC)
RBasin_inflow_mC = [a*b for a,b in zip(RBasin_inflowC,conv_cfs_cms)]
DBasin_inflow_mC = [a*b for a,b in zip(DBasin_inflowC,conv_cfs_cms)]
Wetland_inflow_mC = [a*b for a,b in zip(Wetland_inflowC,conv_cfs_cms)]
Channel_flow_mC = [a*b for a,b in zip(Channel_flowC,conv_cfs_cms)]

# Convert outflow rate from cfs to m3/s
conv_cfs_cms = [0.02832]*len(RBasin_inflowC)
RBasin_outflow_mC = [a*b for a,b in zip(RBasin_outflowC,conv_cfs_cms)]
DBasin_outflow_mC = [a*b for a,b in zip(DBasin_outflowC,conv_cfs_cms)]
Wetland_outflow_mC = [a*b for a,b in zip(Wetland_outflowC,conv_cfs_cms)]

# Convert flooding rate from cfs to m3/s
#conv_cfs_cms = [0.02832]*len(RBasin_inflowC)
#RBasin_flooding_mC = [a*b for a,b in zip(RBasin_floodingC,conv_cfs_cms)]
#Wetland_flooding_mC = [a*b for a,b in zip(Wetland_floodingC,conv_cfs_cms)]
#Wtlnd_bypass_mC = [a*b for a,b in zip(Wtlnd_bp_inflowsC,conv_cfs_cms)]

# Convert depth from ft to m
conv_ft_m = [0.3048]*len(RBasin_inflowC)
RBasin_depth_mC = [a*b for a,b in zip(RBasin_depthC,conv_ft_m)]
DBasin_depth_mC = [a*b for a,b in zip(DBasin_depthC,conv_ft_m)]
Wetland_depth_mC = [a*b for a,b in zip(Wetland_depthC,conv_ft_m)]
Channel_depth_mC = [a*b for a,b in zip(Channel_depthC,conv_ft_m)]

# Calculate outflow load each timestep
conv_mgs_kgs = [0.000001]*len(RBasin_inflowC)
timestep = [5]*len(RBasin_inflowC)
RBasin_loadC = [a*b*c*d*e for a,b,c,d,e in zip(RBasin_concC,RBasin_outflowC,conv_cfs_cms, conv_mgs_kgs,timestep)]
DBasin_loadC = [a*b*c*d*e for a,b,c,d,e in zip(DBasin_concC,DBasin_outflowC,conv_cfs_cms,conv_mgs_kgs,timestep)]
Wetland_loadC = [a*b*c*d*e for a,b,c,d,e in zip(Wetland_concC,Wetland_outflowC,conv_cfs_cms, conv_mgs_kgs,timestep)]
Channel_loadC = [a*b*c*d*e for a,b,c,d,e in zip(Channel_concC,Channel_flowC,conv_cfs_cms,conv_mgs_kgs,timestep)]

# Calculate cumulative load
RBasin_cumloadC = np.cumsum(RBasin_loadC)
DBasin_cumloadC = np.cumsum(DBasin_loadC)
Wetland_cumloadC = np.cumsum(Wetland_loadC)
Channel_cumloadC = np.cumsum(Channel_loadC)

#----------------------------------------------------------------------#

# Print final load released
print("RBasin:", RBasin_cumload[-1])
print("Doyle Basin:", DBasin_cumload[-1])
print("Wetland:", Wetland_cumload[-1])
print("Channel to Outfall:", Channel_cumload[-1]) 
print("RBasin Controlled:", RBasin_cumloadC[-1])
print("Doyle Basin Controlled:", DBasin_cumloadC[-1])
print("Wetland Controlled:", Wetland_cumloadC[-1])
print("Channel Controlled:", Channel_cumloadC[-1])

#----------------------------------------------------------------------#

# Data for flooding line
RB_flood = [6.0960]*len(RBasin_inflowC)
Wt_bypass = [2.8956]*len(RBasin_inflowC)
Wt_flood  = [2.7432]*len(RBasin_inflowC)

# Plot Result
fig, ax = plt.subplots(7,4, constrained_layout=True)
ax[0,0].plot(RBasin_inflow_m, 'k--', linewidth=2)
ax[0,0].plot(RBasin_inflow_mC, linewidth=2, color='#6CC6D1')
ax[0,0].set_xticks([])
#ax[0,0].set_yticks([0,3,6,9])
#ax[0,0].set_yticklabels(["0","3","6","9"])
#ax[0,0].set_ylim(-0.1,9.05)
#ax[0,0].set_xlim(0,190080)
ax[0,0].set_ylabel("Inflow (m³/s)")

ax[1,0].plot(RBasin_conc, 'k--', linewidth=2)
ax[1,0].plot(RBasin_concC, linewidth=2, color='#6CC6D1')
ax[1,0].set_xticks([])
#ax[1,0].set_yticks([0,20,40])
#ax[1,0].set_yticklabels(["0","20","40"])
#ax[1,0].set_ylim(-0.1,50)
#ax[1,0].set_xlim(0,190080)
ax[1,0].set_ylabel("NO (mg/L)")

ax[2,0].set_ylabel("DO (mg/L)")
ax[2,0].set_xticks([])
#ax[2,0].set_yticks([0,5,10])
#ax[2,0].set_yticklabels(["0","5","10"])
#ax[2,0].set_ylim(-0.1,10.5)
#ax[2,0].set_xlim(0,190080)

ax[3,0].plot(RB_flood, linewidth=2, color='#818282')
ax[3,0].plot(RBasin_depth_m, 'k--', linewidth=2)
ax[3,0].plot(RBasin_depth_mC, linewidth=2, color='#6CC6D1')
ax[3,0].set_xticks([])
#ax[3,0].set_yticks([0,3,6])
#ax[3,0].set_yticklabels(["0","3","6"])
#ax[3,0].set_ylim(-0.1,6.5)
ax[3,0].set_xlim(0,225169)
ax[3,0].set_ylabel("Depth (m)")

ax[4,0].plot(RBasin_valveC, linewidth=2, color='#6CC6D1')
ax[4,0].set_xticks([])
ax[4,0].set_ylabel("Valve % Open")
#ax[4,0].set_yticks([0,0.5,1,1.5, 2])
#ax[4,0].set_yticklabels(["0","0.5","1.0","1.5","2.0"])
#ax[4,0].set_ylim(-0.1,2.1)
ax[4,0].set_xlim(0,225169)

ax[5,0].plot(RBasin_outflow_m, 'k--', linewidth=2)
ax[5,0].plot(RBasin_outflow_mC, linewidth=2, color='#6CC6D1')
ax[5,0].set_xticks([])
#ax[5,0].set_yticks([0,4,8])
#ax[5,0].set_yticklabels(["0","4","8"])
#ax[5,0].set_ylim(-0.1,8.5)
ax[5,0].set_xlim(0,225169)
ax[5,0].set_ylabel("Outflow (m³/s)")

ax[6,0].plot(RBasin_cumload, 'k--', linewidth=2)
ax[6,0].plot(RBasin_cumloadC, linewidth=2, color='#6CC6D1')
#ax[6,0].set_yticks([0,1,2])
#ax[6,0].set_yticklabels(["0","1","2"])
#ax[6,0].set_ylim(-0.1,2.25)
ax[6,0].set_xlim(0,225169)
ax[6,0].set_xticks([0,34560,69120,103680,138240,172800,207360])
ax[6,0].set_xticklabels(["0","2","4","6","8","10","12"])
ax[6,0].set_ylabel("Cum. Load (kg)")
ax[6,0].set_xlabel("Time (days)")

ax[0,1].plot(DBasin_inflow_m, 'k--', linewidth=2)
ax[0,1].plot(DBasin_inflow_mC, linewidth=2, color='#3B4D7A')
ax[0,1].set_xticks([])
#ax[0,1].set_yticks([0,3,6,9])
#ax[0,1].set_yticklabels(["0","3","6","9"])
#ax[0,1].set_ylim(-0.1,9.05)
ax[0,1].set_xlim(0,225169)

ax[1,1].plot(DBasin_conc, 'k--', linewidth=2)
ax[1,1].plot(DBasin_concC, linewidth=2, color='#3B4D7A')
ax[1,1].set_xticks([])
#ax[1,1].set_yticks([0,20,40])
#ax[1,1].set_yticklabels(["0","20","40"])
#ax[1,1].set_ylim(-0.1,50)
ax[1,1].set_xlim(0,225169)

ax[2,1].set_xticks([])
#ax[2,1].set_yticks([0,5,10])
#ax[2,1].set_yticklabels(["0","5","10"])
#ax[2,1].set_ylim(-0.1,10.5)
ax[2,1].set_xlim(0,225169)

ax[3,1].plot(Wt_bypass, linewidth=2, color='#818282')
ax[3,1].plot(DBasin_depth_m, 'k--', linewidth=2)
ax[3,1].plot(DBasin_depth_mC, linewidth=2, color='#3B4D7A')
ax[3,1].set_xticks([])
#ax[3,1].set_yticks([0,1.5,3])
#ax[3,1].set_yticklabels(["0","1.5","3"])
#ax[3,1].set_ylim(-0.1,3.25)
ax[3,1].set_xlim(0,225169)

ax[4,1].set_xticks([])
#ax[4,1].set_yticks([0,0.5,1])
#ax[4,1].set_yticklabels(["0","0.5","1.0"])
#ax[4,1].set_ylim(-0.1,1.05)
ax[4,1].set_xlim(0,225169)

ax[5,1].plot(DBasin_outflow_m,'k--', linewidth=2)
ax[5,1].plot(DBasin_outflow_mC, linewidth=2, color='#3B4D7A')
ax[5,1].set_xticks([])
#ax[5,1].set_yticks([0,4,8])
#ax[5,1].set_yticklabels(["0","4","8"])
#ax[5,1].set_ylim(-0.1,8.5)
ax[5,1].set_xlim(0,225169)

ax[6,1].plot(DBasin_cumload, 'k--', linewidth=2)
ax[6,1].plot(DBasin_cumloadC, linewidth=2, color='#3B4D7A')
#ax[6,1].set_yticks([0,1,2])
#ax[6,1].set_yticklabels(["0","1","2"])
#ax[6,1].set_ylim(-0.1,2.25)
ax[6,1].set_xlim(0,225169)
ax[6,1].set_xticks([0,34560,69120,103680,138240,172800,207360])
ax[6,1].set_xticklabels(["0","2","4","6","8","10","12"])
ax[6,1].set_xlabel("Time (days)")

ax[0,2].plot(Wetland_inflow_m, 'k--', linewidth=2)
ax[0,2].plot(Wetland_inflow_mC, linewidth=2, color='#B08CA1')
ax[0,2].set_xticks([])
#ax[0,2].set_yticks([0,3,6,9])
#ax[0,2].set_yticklabels(["0","3","6","9"])
#ax[0,2].set_ylim(-0.1,9.05)
ax[0,2].set_xlim(0,225169)

ax[1,2].plot(Wetland_conc, 'k--', linewidth=2)
ax[1,2].plot(Wetland_concC, linewidth=2, color='#B08CA1')
ax[1,2].set_xticks([])
#ax[1,2].set_yticks([0,20,40,60,80])
#ax[1,2].set_yticklabels(["0","20","40","60","80"])
#ax[1,2].set_ylim(-0.1,100)
ax[1,2].set_xlim(0,225169)

ax[2,2].plot(Wetland_DO, 'k--', linewidth=2)
ax[2,2].plot(Wetland_DOC, linewidth=2, color="#B08CA1")
ax[2,2].set_xticks([])
#ax[2,2].set_yticks([0,5,10,15,20])
#ax[2,2].set_yticklabels(["0","5","10","15","20"])
#ax[2,2].set_ylim(-0.1,25)
ax[2,2].set_xlim(0,225169)

ax[3,2].plot(Wt_flood, linewidth=2, color='#818282')
ax[3,2].plot(Wetland_depth_m, 'k--', linewidth=2)
ax[3,2].plot(Wetland_depth_mC, linewidth=2, color='#B08CA1')
ax[3,2].set_xticks([])
#ax[3,2].set_yticks([0,1.5,3])
#ax[3,2].set_yticklabels(["0","1.5","3"])
#ax[3,2].set_ylim(-0.1,3.25)
ax[3,2].set_xlim(0,225169)

ax[4,2].plot(Wetland_valveC, linewidth=2, color='#B08CA1')
ax[4,2].set_xticks([])
ax[4,2].set_yticks([0,0.5,1,1.5,2])
#ax[4,2].set_yticklabels(["0","0.5","1.0","1.5","2"])
#ax[4,2].set_ylim(-0.1,2.1)
ax[4,2].set_xlim(0,225169)

ax[5,2].plot(Wetland_outflow_m, 'k--', linewidth=2)
ax[5,2].plot(Wetland_outflow_mC, linewidth=2, color='#B08CA1')
ax[5,2].set_xticks([])
#ax[5,2].set_yticks([0,4,8])
#ax[5,2].set_yticklabels(["0","4","8"])
#ax[5,2].set_ylim(-0.1,8.5)
ax[5,2].set_xlim(0,225169)

ax[6,2].plot(Wetland_cumload, 'k--', linewidth=2)
ax[6,2].plot(Wetland_cumloadC, linewidth=2, color='#B08CA1')
#ax[6,2].set_yticks([0,1,2])
#ax[6,2].set_yticklabels(["0","1","2"])
#ax[6,2].set_ylim(-0.1,2.25)
ax[6,2].set_xlim(0,225169)
ax[6,2].set_xticks([0,34560,69120,103680,138240,172800,207360])
ax[6,2].set_xticklabels(["0","2","4","6","8","10","12"])
ax[6,2].set_xlabel("Time (days)")

ax[0,3].plot(Channel_flow_m, 'k--', linewidth=2)
ax[0,3].plot(Channel_flow_mC, linewidth=2, color='#695580')
ax[0,3].set_xticks([])
#ax[0,3].set_yticks([0,3,6,9])
#ax[0,3].set_yticklabels(["0","3","6","9"])
#ax[0,3].set_ylim(-0.1,9.05)
ax[0,3].set_xlim(0,225169)

ax[1,3].plot(Channel_conc, 'k--', linewidth=2)
ax[1,3].plot(Channel_concC, linewidth=2, color='#695580')
ax[1,3].set_xticks([])
#ax[1,3].set_yticks([0,20,40])
#ax[1,3].set_yticklabels(["0","20","40"])
#ax[1,3].set_ylim(-0.1,50)
ax[1,3].set_xlim(0,225169)

ax[2,3].set_xticks([])
#ax[2,3].set_yticks([0,5,10])
#ax[2,3].set_yticklabels(["0","5","10"])
#ax[2,3].set_ylim(-0.1,10.5)
ax[2,3].set_xlim(0,225169)

ax[3,3].plot(Channel_depth_m, 'k--', linewidth=2)
ax[3,3].plot(Channel_depth_mC, linewidth=2, color='#695580')
ax[3,3].set_xticks([])
#ax[3,3].set_yticks([0,1.5,3])
#ax[3,3].set_yticklabels(["0","1.5","3"])
#ax[3,3].set_ylim(-0.1,3.25)
ax[3,3].set_xlim(0,225169)

ax[4,3].set_xticks([])
#ax[4,3].set_yticks([0,0.5,1])
#ax[4,3].set_yticklabels(["0","0.5","1.0"])
#ax[4,3].set_ylim(-0.1,1.05)
ax[4,3].set_xlim(0,225169)

ax[5,3].plot(Channel_flow_m, 'k--', linewidth=2)
ax[5,3].plot(Channel_flow_mC, linewidth=2, color='#695580')
ax[5,3].set_xticks([])
#ax[5,3].set_yticks([0,4,8])
#ax[5,3].set_yticklabels(["0","4","8"])
#ax[5,3].set_ylim(-0.1,8.5)
ax[5,3].set_xlim(0,225169)

ax[6,3].plot(Channel_cumload, 'k--', linewidth=2)
ax[6,3].plot(Channel_cumloadC, linewidth=2, color='#695580')
#ax[6,3].set_yticks([0,1,2])
#ax[6,3].set_yticklabels(["0","1","2"])
#ax[6,3].set_ylim(-0.1,2.25)
ax[6,3].set_xlim(0,225169)
ax[6,3].set_xticks([0,34560,69120,103680,138240,172800,207360])
ax[6,3].set_xticklabels(["0","2","4","6","8","10","12"])
ax[6,3].set_xlabel("Time (days)")
plt.savefig('NOresults.eps')
plt.show()
