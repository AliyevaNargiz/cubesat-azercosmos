# Cubesat GEO Mission Planner

Plan and analyze a cubesat transfer from Kourou (French Guiana) to a GEO slot between HORIZONS‑1 and INMARSAT‑3 F3 using Orekit 13.1 and Hipparchus 4.0.1.

## Contents
- `GeoMissionPlanner.java`: End‑to‑end GEO mission planning (neighbors, GTO, apogee burn ΔV, drift, GMAT‑ready output)
- `LaunchAnalysis.java`: Simpler launch and transfer analysis with summary outputs
- `run_mission.bat`, `run_launch.bat`, `run.bat`: Windows helpers to run the programs
- `*.jar`: Required Orekit and Hipparchus libraries

## Prerequisites
- Java JDK 11+ (`java -version`)
- Orekit data bundle (EOP, time scales, etc.)
  - Download the official dataset zip and unzip it, e.g. into `orekit-data-main/`.
  - Source: [Orekit data repository](https://gitlab.orekit.org/orekit/orekit-data)

## Quick Start (Windows / PowerShell)
1) Open this folder in PowerShell:
```powershell
cd C:\Users\hp\Desktop\Cubesat
```

2) Ensure Orekit data is present (folder like `orekit-data-main`):
```powershell
ls orekit-data-main
```

3) Compile (already compiled if `*.class` exist, but safe to run):
```powershell
javac -cp "hipparchus-core-4.0.1.jar;hipparchus-geometry-4.0.1.jar;hipparchus-ode-4.0.1.jar;hipparchus-optim-4.0.1.jar;hipparchus-stat-4.0.1.jar;hipparchus-fft-4.0.1.jar;hipparchus-filtering-4.0.1.jar;hipparchus-fitting-4.0.1.jar;hipparchus-samples-4.0.1.jar;hipparchus-clustering-4.0.1.jar;orekit-13.1.jar" GeoMissionPlanner.java LaunchAnalysis.java
```

4) Run the full mission planner:
```powershell
./run_mission.bat orekit-data-main
```

5) (Optional) Run the simpler launch analysis:
```powershell
./run_launch.bat orekit-data-main
```

If you prefer setting an environment variable instead of passing the path each time:
```powershell
setx OREKIT_DATA_PATH "C:\Users\hp\Desktop\Cubesat\orekit-data-main"
# then you can run without the argument (GeoMissionPlanner still accepts it)
```

## What You Get
- Neighbor satellite longitudes (HORIZONS‑1 and INMARSAT‑3 F3) near your chosen launch epoch
- GTO design (a, e, i; perigee/apogee altitudes)
- Time to apogee and apogee burn ΔV (circularization + plane change)
- GEO drift plan to target longitude (Δa, drift rate, time, small ΔV set/stop)
- Safety separation check vs. neighbors at arrival
- GMAT‑ready final GEO orbital elements

## Change Inputs
Adjust these in `GeoMissionPlanner.java`:
- Launch site (Kourou defaults provided)
- Launch epoch (UTC)
- Initial GTO perigee/apogee and inclination
- Target GEO longitude (slot center)
- Final GEO inclination target
- Whether to include perigee injection from LEO (`includePerigeeInjection`)

## Troubleshooting
- Error: `no IERS UTC-TAI history data loaded`
  - Ensure the Orekit data folder exists and is passed as an argument, or set `OREKIT_DATA_PATH`.
  - The folder should contain files like `tai-utc.dat`, EOP files under `Earth-Orientation-Parameters/`, etc.
- PowerShell quoting or wildcard issues
  - Use the provided `run_*.bat` scripts to avoid classpath quoting problems.
- Garbled degree symbols in console
  - This is a console encoding issue only; computations are unaffected.

## Repository
GitHub: [`cubesat-azercosmos`](https://github.com/AliyevaNargiz/cubesat-azercosmos)

## Notes
- Large assets are intentionally not tracked (see `.gitignore`). Keep `*.jar` locally or switch to a dependency manager (e.g., Maven/Gradle) if you prefer reproducible builds.
- This planner uses simplified assumptions for clarity (e.g., Keplerian propagation for some steps). For flight‑ready analysis, include higher‑fidelity modeling and operations constraints.
