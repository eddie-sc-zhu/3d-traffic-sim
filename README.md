# 3d-traffic-sim
A 3 Dimensional Dynamic Traffic Assignment Engine built upon the Unity Game Engine


# How To Access Source Code
1) Download Unity 2019.4.3f1 from here: https://unity3d.com/get-unity/download
2) Download C# .NET 4.x or higher (due to use of dynamic parameter)
3) Clone this project to your desired file location
4) Done
# Updates
Currently Work in progress:
- [x] GMNS data parser
- [x] DTA-lite functionality
- [x] User-equilibrium traffic assignment (Frank-Wolfe / gradient projection)
- [ ] 3D Visualization
- [ ] Equity-based screening/statistics
- [ ] Blender-osm/simviz narrow simulation

# Traffic Assignment
`Assets/Data/DLSim.py` solves the static traffic assignment problem to user
equilibrium (Wardrop's first principle, Beckmann formulation with BPR volume
delay functions) before seeding the mesoscopic simulation, so agents drive
equilibrium routes instead of naive shortest paths. Configure via the
constants at the top of `DLSim.py`:
- `ASSIGNMENT_MODE`: `'UE_FW'` (Frank-Wolfe, link-based, default),
  `'UE_GP'` (gradient projection, path-based), or `'AON'` (legacy one-shot
  all-or-nothing shortest paths)
- `UE_MAX_ITERATIONS`, `UE_RELATIVE_GAP_TOLERANCE`: convergence controls
- `UE_DEMAND_MULTIPLIER`: scales OD demand to stress the network

The solver lives in `Assets/Data/traffic_assignment.py` and writes
`ue_link_performance.csv` (equilibrium volumes, v/c ratios, congested travel
times), `ue_route_assignment.csv` (equilibrium path flows per OD pair) and
`ue_convergence.csv` (relative gap per iteration).

# How To Build (Not Implemented)
Building will be in a single executable file. In the source code there will be a folder called "builds." Click that and open the executable file. It should open the program and be able to run.
 
# References

Sebastian Lague Path-Creator: https://github.com/SebLague/Path-Creator

