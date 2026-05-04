# UR5 Robotic Arm Dynamics, Control, and Environmental Interaction

MATLAB simulation environment for studying the dynamics, control, trajectory tracking, impedance control, and environmental interaction of a UR5 robotic manipulator.

<img width="1562" height="988" alt="ur5_pick_place_dynamics_control" src="https://github.com/user-attachments/assets/d406826a-8fd1-4c0a-98e5-e1bfd4dc5edb" />


## Overview

This repository provides MATLAB scripts for simulating and comparing multiple robotic-arm control strategies for a 6-DOF UR5 manipulator. The project includes generated robot dynamics functions, parameter definitions, controller simulations, trajectory tracking experiments, impedance-control scenarios, inverse-kinematics task animation, full nonlinear dynamics simulation, and constrained environmental interaction tasks.

The main focus areas are:

- Rigid-body dynamics and kinematics of the UR5 manipulator
- Joint-space regulation and tracking control
- PD, PID, gravity-compensated, and feedback-linearization-based controllers
- Numerical inverse kinematics for task-space motion generation
- Full nonlinear dynamic simulation with computed-torque control
- External payload/load compensation during pick-and-place motion
- Parameter uncertainty studies
- Cartesian impedance control under external forces
- Hybrid force/motion control during contact with a plane
- Constrained path tracking, including sphere-constrained motion
- Reduced dynamic parameterization analysis

## Repository Structure

```text
.
│
├── maple_gen/
│   ├── UR5_C.m
│   ├── UR5_G.m
│   ├── UR5_M.m
│   ├── UR5_MCG.m
│   ├── UR5_fdyn.m
│   ├── UR5_fkall.m
│   ├── UR5_fkine.m
│   ├── UR5_h.m
│   ├── UR5_idyn.m
│   ├── UR5_jacobian_geometric.m
│   └── symbolic_robot_dyns.mw
│
├── utils/
│   ├── figureoptscall.m
│   ├── inverse_kinematics.m
│   └── saveFigureAsPDF.m
│
├── UR5_params.m
├── rigidbodytree_test.m
├── test_gen_codes.m
├── run_reduced_parameterization.m
├── run_latex_reduced_parameterization.m
├── run_ur5_pd_comparison.m
├── run_ur5_pid_sweep.m
├── run_ur5_pd_iterative_gravity.m
├── run_ur5_fb_tracking_known.m
├── run_ur5_fb_tracking_uncertain.m
├── run_ur5_cartesian_impedance.m
├── run_ur5_hybrid_plane_force_motion.m
├── run_ur5_sphere_constrained_tracking.m
├── run_ur5_pick_place_ik.m
└── run_ur5_pick_place_dynamics_control.m
```

## Main Components

### `UR5_params.m`

Defines the robot geometric and dynamic parameters and packs them into the parameter vector `Pi` used by the generated dynamics functions.

The parameter set includes:

- Link geometry
- Link masses
- Center-of-mass positions
- Inertia tensor terms
- Gravity

### `maple_gen/`

Contains the generated MATLAB functions for UR5 dynamics and kinematics, together with the original Maple worksheet used for symbolic model development.

The file `symbolic_robot_dyns.mw` is the Maple worksheet containing the symbolic derivation of the robot dynamics. It serves as the source worksheet for deriving the analytical model and generating the MATLAB functions used by the simulations.

The generated MATLAB functions include:

- Mass matrix
- Coriolis/centrifugal terms
- Gravity vector
- Forward dynamics
- Inverse dynamics
- Forward kinematics
- Full link transforms
- Geometric Jacobian

These functions are used by the simulation scripts and should remain on the MATLAB path.

### `utils/`

Contains helper functions for numerical inverse kinematics, figure formatting, and exporting plots to PDF.

The main utility files are:

- `inverse_kinematics.m`: reusable damped least-squares numerical inverse-kinematics solver
- `figureoptscall.m`: common figure formatting and LaTeX-style plot configuration
- `saveFigureAsPDF.m`: PDF export utility used by the simulation scripts

### `assets/`

Contains media files used by the README, including the exported pick-and-place simulation video.

## Simulation Scripts

| Script | Purpose |
|---|---|
| `test_gen_codes.m` | Tests the generated dynamics and kinematics functions. |
| `rigidbodytree_test.m` | Compares or validates the generated model against a MATLAB rigid-body-tree style model. |
| `run_reduced_parameterization.m` | Performs reduced dynamic parameterization analysis. |
| `run_latex_reduced_parameterization.m` | Generates LaTeX-oriented reduced-parameterization output. |
| `run_ur5_pd_comparison.m` | Compares PD, gravity-compensated PD, approximate-gravity PD, desired-gravity PD, and PID regulation controllers. |
| `run_ur5_pid_sweep.m` | Sweeps integral gain values for PID regulation and compares tracking/error metrics. |
| `run_ur5_pd_iterative_gravity.m` | Studies PD control with iterative gravity-related compensation. |
| `run_ur5_fb_tracking_known.m` | Compares feedforward, variable-PD, feedback-linearization, and PID feedback-linearization tracking under known model parameters. |
| `run_ur5_fb_tracking_uncertain.m` | Studies tracking performance under model uncertainty. |
| `run_ur5_cartesian_impedance.m` | Simulates Cartesian impedance control with and without external force disturbances. |
| `run_ur5_hybrid_plane_force_motion.m` | Simulates hybrid force/motion control while interacting with a rigid plane. |
| `run_ur5_sphere_constrained_tracking.m` | Simulates constrained tracking on a spherical surface. |
| `run_ur5_pick_place_ik.m` | Generates and animates a pick-and-place task using numerical inverse kinematics. |
| `run_ur5_pick_place_dynamics_control.m` | Simulates the pick-and-place task using full nonlinear dynamics, computed-torque control, and payload torque compensation. |

## Requirements

The project is written in MATLAB.

Recommended tools and toolboxes:

- MATLAB base environment
- Maple base environment, if regenerating or modifying symbolic dynamics
- Robotics System Toolbox, if using or modifying rigid-body-tree validation scripts

No external Python or C++ dependencies are required for the provided MATLAB simulations.

## Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/Hosnooo/UR50-Robotic-Arm-Dynamics-Control-Enviromental-Interaction.git
cd UR50-Robotic-Arm-Dynamics-Control-Enviromental-Interaction
```

### 2. Open MATLAB in the repository folder

Start MATLAB and set the current folder to the repository root.

### 3. Add required folders to the path

Most scripts already include:

```matlab
addpath('maple_gen')
addpath('utils')
```

You can also add them manually:

```matlab
addpath('maple_gen');
addpath('utils');
```

### 4. Run a test script

To verify that the generated model functions are accessible:

```matlab
test_gen_codes
```

## Example Usage

### Run a PD/PID regulation comparison

```matlab
run_ur5_pd_comparison
```

This script compares several joint-space regulation controllers and saves simulation results, metrics, and plots.

### Run a PID integral-gain sweep

```matlab
run_ur5_pid_sweep
```

This script evaluates how different integral gain scales affect regulation error, settling time, torque demand, and control energy.

### Run known-parameter trajectory tracking

```matlab
run_ur5_fb_tracking_known
```

This script compares several tracking controllers, including feedforward PD, variable PD, feedback linearization with PD/feedforward, and feedback linearization with PID/feedforward.

### Run Cartesian impedance control

```matlab
run_ur5_cartesian_impedance
```

This script simulates Cartesian impedance behavior for nominal, soft, and stiff impedance settings under external force disturbances.

### Run hybrid force/motion plane interaction

```matlab
run_ur5_hybrid_plane_force_motion
```

This script simulates motion control in the tangent directions of a plane while regulating the contact force in the plane-normal direction.

### Run inverse-kinematics pick-and-place animation

```matlab
run_ur5_pick_place_ik
```

This script generates a smooth pick-and-place task in Cartesian space, solves the corresponding joint trajectory using the numerical inverse-kinematics utility, animates the UR5 motion, and saves a video.

### Run full-dynamics pick-and-place control

```matlab
run_ur5_pick_place_dynamics_control
```

This script first generates a pick-and-place reference trajectory using numerical inverse kinematics. It then simulates the full nonlinear UR5 dynamics under a computed-torque tracking controller. During the carrying phase, an external payload torque is included and compensated through the controller.

## Pick-and-Place Study

The pick-and-place task is split into two levels:

### Kinematic level

The script `run_ur5_pick_place_ik.m` uses numerical inverse kinematics to answer:

```text
What joint trajectory allows the end-effector to follow the desired pick-and-place path?
```

It uses:

- `UR5_fkine.m`
- `UR5_fkall.m`
- `UR5_jacobian_geometric.m`
- `utils/inverse_kinematics.m`

### Dynamic level

The script `run_ur5_pick_place_dynamics_control.m` uses the IK-generated joint trajectory as the reference for a full nonlinear dynamic simulation.

The controller is a computed-torque tracking controller of the form:

```text
tau = M(q) v + h(q,dq) + G(q)
```

where:

```text
v = ddq_d + Kd(dq_d - dq) + Kp(q_d - q)
```

When the payload is carried, the payload force is mapped into joint space as:

```text
tau_ext = J_v(q)^T F_ext
```

and compensated in the motor command using:

```text
tau_motor = tau - tau_ext
```

The plant then receives:

```text
tau_total = tau_motor + tau_ext
```

which keeps the simulation physically meaningful while allowing compensation of the known external payload load.

## Output Files

Each simulation creates a dedicated results folder. Depending on the script, outputs may include:

- CSV time-series logs
- Summary tables
- MAT files containing simulation structures and metrics
- MATLAB figures
- PDF figures, when the figure export utilities are available
- MP4 or AVI videos for animation scripts

Example result folders include:

```text
ur5_pd_comparison_results/
ur5_pid_sweep_results/
ur5_fb_tracking_results/
ur5_fb_uncertainty_results/
ur5_cartesian_impedance_results/
ur5_hybrid_plane_results/
ur5_pick_place_dynamics_results/
```

## Typical Metrics

The scripts compute and save metrics such as:

- Final position error
- RMS position error
- Peak position error
- Final velocity error
- RMS velocity error
- Peak torque norm
- RMS torque norm
- Control energy
- Settling time
- Contact-force error
- Constraint violation magnitude
- Joint tracking error
- End-effector tracking error

## Control Methods Included

The repository includes simulations for several control approaches.

### Joint-Space Regulation

- PD control
- PD with exact gravity compensation
- PD with approximate gravity compensation
- PD with desired-configuration gravity compensation
- PID control

### Joint-Space Tracking

- Feedforward plus PD
- Feedforward plus variable PD
- Feedback linearization with PD and feedforward
- Feedback linearization with PID and feedforward
- Computed-torque tracking control

### Cartesian and Interaction Control

- Numerical inverse kinematics
- Cartesian impedance control
- External force disturbance response
- External payload torque compensation
- Hybrid force/motion control
- Plane-constrained contact simulation
- Sphere-constrained trajectory tracking

## Notes on Symbolic and Generated Dynamics

The symbolic robot model is developed in the Maple worksheet `maple_gen/symbolic_robot_dyns.mw`. This worksheet documents the symbolic derivation process and is used to generate the analytical MATLAB dynamics and kinematics functions.

The MATLAB functions in `maple_gen/` are generated UR5 model functions and are expected to be called by the simulation scripts. Avoid editing the generated `.m` files manually unless you are intentionally changing the robot model or regenerating the symbolic dynamics from the Maple worksheet.

The parameter vector `Pi` created by `UR5_params.m` must match the parameter ordering expected by the generated functions.

## Recommended Workflow

1. Run `test_gen_codes.m` to verify generated functions.
2. Run `run_ur5_pd_comparison.m` to test basic regulation controllers.
3. Run `run_ur5_pid_sweep.m` to tune integral gains.
4. Run `run_ur5_fb_tracking_known.m` for full-model tracking comparison.
5. Run `run_ur5_fb_tracking_uncertain.m` to study robustness.
6. Run `run_ur5_pick_place_ik.m` to visualize a task-space pick-and-place trajectory.
7. Run `run_ur5_pick_place_dynamics_control.m` to simulate full nonlinear dynamics with computed-torque control.
8. Run Cartesian/environment-interaction scripts for advanced simulations.

## Citation

If this repository is used for academic work, cite the repository and any related project report, thesis, or paper associated with the model and controller design.

```bibtex
@software{ur5_robotic_arm_dynamics_control_environmental_interaction,
  title  = {UR5 Robotic Arm Dynamics, Control, and Environmental Interaction},
  author = {Mohssen Elshaar and Nathan Lablanc},
  year   = {2026},
  url    = {https://github.com/Hosnooo/UR50-Robotic-Arm-Dynamics-Control-Enviromental-Interaction}
}
```
