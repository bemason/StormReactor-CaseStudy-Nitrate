# Nitrate SWMM Simulations
Python scripts for running nitrate simulations for StormReactor research article

## Files
1. ΝΟ.inp: SWMM input file for nitrate simulations based on real-world inspired stormwater network.
2. NO_check.inp: Same as #1 but shorter simulation time (3 days vs 12 days). 
3. toolbox_NO.py: Primary Python script used for all nitrate methods, both controlled and uncontrolled scenarios, run for StormReactor research paper.
4. toolbox_NO_timer_RTC.py: Python script to check time requirement for running the nitrate simulation with real-time control.
5. toolbox_NO_timer_noRTC.py: Python script to check time requirement for running the nitrate simulation without real-time control.
6. toolbox_NO_timer_noWQ.py: Python script to check time requirement for running the simulation without the nitrate modele or real-time control.
7. NO_check.py: Python script to compare steady state analytical nitrate concentration to StormReactor's nitrate model.

## License
GNU General Public License v3.0

## Status
These Python scripts use wq_toolbox, the prerequisite package to StormReactor. These files will be updated once StormReactor is finalized.
