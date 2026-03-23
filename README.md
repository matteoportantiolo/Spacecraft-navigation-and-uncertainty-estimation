# Spacecraft Navigation and Uncertainty Estimation

## Overview
This project focuses on uncertainty propagation and state estimation techniques applied to astrodynamics problems, with applications to Earth-Moon transfers, Earth-orbiting satellites, and lunar navigation.

The work combines linear and nonlinear estimation methods, including covariance propagation, Monte Carlo simulations, batch least-squares estimation, and sequential filtering using Unscented Kalman Filters (UKF).

## Mission Context
The project addresses three main navigation scenarios:
- Uncertainty propagation in nonlinear multi-body dynamics (Earth-Moon system)
- Orbit determination of an Earth observation satellite (SMOS) using ground stations
- Navigation and localization in a lunar environment using inter-satellite measurements

Key objectives include:
- Accurate uncertainty propagation in nonlinear dynamics
- Optimal use of measurement networks for orbit determination
- Robust sequential estimation for real-time navigation

## System Architecture

### Dynamical Models
Different dynamical models are considered depending on the scenario:
- Planar Bicircular Restricted Four-Body Problem (PBRFBP)
- Keplerian two-body motion
- J2-perturbed orbital dynamics
- SGP4 propagation for realistic LEO trajectories

### Measurement Models
- Ground station measurements: Azimuth, Elevation, Range
- Inter-satellite measurements: Relative range
- Direct position measurements (lunar navigation service)

Noise is modeled as Gaussian with known covariance.

## Methods

### Uncertainty Propagation
Three approaches are implemented and compared:
- Linearized covariance propagation (STM-based)
- Unscented Transform (UT)
- Monte Carlo simulation

The comparison highlights the limitations of linear methods in strongly nonlinear dynamics

### Batch Estimation
A nonlinear least-squares approach is used for orbit determination:
- Levenberg-Marquardt algorithm (`lsqnonlin`)
- Weighted residual formulation
- Comparison of:
  - Single station vs multi-station tracking
  - Keplerian vs J2 dynamics

### Sequential Estimation
An Unscented Kalman Filter (UKF) is implemented for:
- Orbiter state estimation
- Joint estimation of orbiter state and lunar lander coordinates

The UKF propagates sigma points through nonlinear dynamics and measurement models to capture non-Gaussian effects.

### Trade-off Analysis
Ground station selection is optimized based on:
- Estimation accuracy (semi-major axis and inclination uncertainty)
- Operational cost constraints

## Results

### Uncertainty Propagation
- UT provides results close to Monte Carlo simulations
- Linearized methods underestimate uncertainty in nonlinear regimes

### Orbit Determination
- Multi-station tracking significantly improves accuracy
- Including J2 perturbation reduces estimation error by orders of magnitude

### Trade-off Analysis
- Best cost-performance trade-off achieved with Kourou–Svalbard combination
- Long-term operations favor polar stations due to orbital geometry

### Sequential Filtering
- UKF successfully reduces estimation error over time
- Errors remain within 3σ bounds, validating filter consistency
- Joint estimation of orbiter and lander is feasible with limited measurements

## Implementation
The project includes:
- Numerical propagation of orbital dynamics
- STM integration for linear covariance propagation
- Monte Carlo simulation framework
- Nonlinear least-squares solver for orbit determination
- Unscented Kalman Filter implementation
- Trade-off and mission analysis tools

## Key Concepts
- Orbit determination
- Uncertainty propagation
- Unscented Transform (UT)
- Unscented Kalman Filter (UKF)
- Monte Carlo methods
- Nonlinear least squares
- Ground station visibility analysis
- J2 perturbation effects

## Author
Matteo Portantiolo  
MSc Space Engineering – GNC
