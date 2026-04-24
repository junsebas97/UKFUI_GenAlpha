# Unscented Kalman Filter with Uncertain Inputs based on the Generalised- $\alpha$ method (UKFUI-G $\alpha$)
This repository implements the UKFUI-G $\alpha$ [1]. This filter estimates states and uncertain inputs of nonlinear structures using Gaussian random walks, the unscented transform, and the generalised- $\alpha$ method with $\rho_\infty = 1$.

The repository contains three examples on an 8-DOF shear chain with Bouc-Wen springs:

- [a numerical validation and comparison with the UKF-UI [1] and UKF-UL [2]](./numerical_validation.m),
- [a sensitivity analysis on the noise level](./sensitivity_measurement_noise.m), and
- [a sensitivity analysis on the number of measurements](./sensitivity_number_measurement.m).


## Bibliography:
[1] J. S. Delgado Trujillo, J. Tott-Buswell, S. Jalbi, J. Hilton, M. Pandey, L. J. Prendergast, A nonlinear Bayesian filter with uncertain inputs based on the generalised- $\alpha$ method [Forthcoming].

[2] Y. Lei, D. Xia, K. Erazo, S. Nagarajaiah, A novel unscented Kalman filter for recursive state-input-system identification of nonlinear systems, Mechanical Systems and Signal Processing 127 (2019) 120–135. doi:10.1016/j.ymssp.2019.03.013.

[3] J. S. Delgado Trujillo, J. Tott-Buswell, S. Jalbi, J. Hilton, M. Pandey, L. J. Prendergast, A nonlinear Bayesian filter for structural systems subjected to uncertain loads, Journal of Sound and Vibration. [In press]