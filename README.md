# Unscented Kalman Filter with Uncertain Inputs based on the Generalised-$\alpha$ method (UKFUI-G $\alpha$)
This repository implements the UKFUI-G $\alpha$ [CITE PENDING]. This filter estimates states and uncertain inputs of nonlinear structures using Gaussian random walks, the unscented transform, and the generalised-$\alpha$ method with $\rho_\infty = 1$.

The repository contains three examples on an 8-DOF shear chain with Bouc-Wen springs:

- [a numerical validation and comparison with the UKF-UI and UKF](./numerical_validation.m),
- [a sensitivity analysis on the noise level](./sensitivity_measurement_noise.m), and
- [a sensitivity analysis on the number of measurements](./sensitivity_number_measurement.m).
