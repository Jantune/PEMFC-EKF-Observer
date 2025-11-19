# Cell-by-Cell Temperature Observer for Commercial PEMFC Stacks

## Overview
This repository contains the **open-source MATLAB implementation** for the paper:

**"Cell-by-Cell Temperature Observer for Commercial PEMFC Stacks: A Thermal-Hydraulic Approach"**

The project implements a physics-based Extended Kalman Filter (EKF) designed to estimate the internal temperature distribution of a large-scale PEMFC stack in real-time, using only standard onboard sensors.

## Key Features
- ✅ **Analytical Tridiagonal Jacobian**: Reduces computational complexity to $O(n^2)$
- ✅ **Thermal-Hydraulic Coupling**: Integrated flow distribution and 1D thermal model
- ✅ **Fully Data-Adaptive**: 
  - Automatically detects cell count ($N$) from input data
  - Automatically calculates time step ($\Delta t$) from time vector
  - No hardcoded dimensions or snapshot times
- ✅ **Centralized Configuration**: All physical constants in `LoadStackConfig.m`
- ✅ **Robustness**: Validated against coolant interruption faults and dynamic loads
- ✅ **Benchmarking**: Includes ML baseline comparisons (Ridge, GPR)

## Project Structure

```
PEMFC-EKF-Observer/
├── Run_All.m                  # Main entry point
├── Data/
│   ├── LoadData.m             # Data loading interface (populate with your data)
│   ├── LoadStackConfig.m      # Centralized physical constants
│   └── GenerateSampleData.m   # Synthetic test data generator
├── Models/
│   ├── HydraulicModel.m       # Coolant flow distribution
│   └── ThermalModel.m         # 1D unsteady thermal physics
├── EKF/
│   ├── ExtendedKalmanFilter.m # Core observer with analytical Jacobian
│   └── Run_EKF_Estimation.m   # Main estimation script
├── Optimization/
│   ├── SimplePSO.m            # Particle Swarm Optimization
│   ├── Identify_Hydraulic_Parameters.m
│   ├── Identify_Thermal_Parameters.m
│   └── Optimize_QR_Matrices.m
├── Analysis/
│   ├── TimePerformanceAnalysis.m  # Computational benchmarks
│   └── RobustnessAnalysis.m       # Parameter sensitivity
└── ML_Baselines/
    ├── TrainResidualBaseline.m
    └── Compare_ML_Models.m
```

## Input Data Requirements

Modify `Data/LoadData.m` to load your experimental data. The function must return a MATLAB `struct` with:

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `time` | $T \times 1$ | s | Time vector |
| `I_st` | $T \times 1$ | A | Stack current |
| `flow_rate` | $T \times 1$ | L/min | Total coolant flow rate |
| `T_in` | $T \times 1$ | °C | Coolant inlet temperature |
| `T_out` | $T \times 1$ | °C | Coolant outlet temperature (measurement for EKF) |
| `V_cell` | $T \times N$ | V | Individual cell voltages |
| `P_in` | $T \times 1$ | kPa | Inlet pressure (optional, for hydraulic ID) |
| `P_out` | $T \times 1$ | kPa | Outlet pressure (optional, for hydraulic ID) |
| `T_dist` | $T \times N$ | °C | True cell temperatures (optional, for validation) |

**Important Notes:**
- ✅ **$N$ is auto-detected** from `size(V_cell, 2)` - no hardcoding required
- ✅ **$\Delta t$ is auto-calculated** from `mean(diff(time))` 
- ✅ **Time range is flexible** - any duration and sampling rate supported
- ✅ **Snapshot times are adaptive** - ML training uses intelligent fallback if paper-specific times (635s, 1349s, 1708s) not found

## Physical Constants Configuration

All stack-specific constants are centralized in:
```matlab
Data/LoadStackConfig.m
```
This includes:
- Hydraulic geometry (channel dimensions, resistances)
- Material properties (viscosity, density)
- Thermal/electrochemical constants (E₀, Tenv, Aed)
- Default identified parameters (C, R, h, α, β)
- Default EKF tuning (Q, R matrices)

**To adapt to your stack:** Modify `LoadStackConfig.m` only - no other files need editing.

## Requirements
- MATLAB R2020b or newer
- Optimization Toolbox (recommended, though `SimplePSO.m` included as fallback)
- Statistics and Machine Learning Toolbox (for GPR baseline only)

## Usage

### Quick Start (Synthetic Data)
```matlab
cd PEMFC-EKF-Observer
Run_All  % Runs full pipeline with sample data
```

### With Real Data
1. Edit `Data/LoadData.m` to load your experimental data
2. (Optional) Update `Data/LoadStackConfig.m` for your stack geometry
3. Run `Run_All.m`

### Individual Modules
```matlab
% Parameter Identification
Identify_Thermal_Parameters()
Identify_Hydraulic_Parameters()

% EKF Estimation
Run_EKF_Estimation()

% Analysis
run('Analysis/TimePerformanceAnalysis.m')
RobustnessAnalysis()

% ML Comparison
Compare_ML_Models()
```

## Adaptivity Examples

**Example 1: Different Stack Size**
- Your data: 200 cells → Code automatically uses $N=200$
- Your data: 50 cells → Code automatically uses $N=50$

**Example 2: Different Time Range**
- Your data: 0-5000s @ 0.5s → Code adapts to 10000 steps @ 0.5s
- Your data: 100-300s @ 2s → Code adapts to 100 steps @ 2s

**Example 3: Different Validation Strategy**
- Has 635s, 1349s, 1708s → Uses paper validation strategy
- No specific times → Automatically uses 70/30 train/val split

## Citation
If you use this code in your research, please cite:
```
[Updating]
```

## License
MIT License

## Contact
For questions or issues, please open an issue on [GitHub](https://github.com/Jantune/PEMFC-EKF-Observer).
