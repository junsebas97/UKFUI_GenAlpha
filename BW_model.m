function [xi_i, r_i, K_BW] = BW_model(xi_im1, v_im1, v_i, u_i, k, alpha, ...
                                       beta, gamma, n, dt)
%{
this function is the discrete Bouc-Wen model. it calculates the new BW
displacement using the 4th-order Runge-Kutta method.

INPUTS:
xi_im1: initial Bouc-Wen displacement   [m]
v_im1:  initial DOFs velocities         [m/s]
v_i:    end DOFs velocities             [m/s]
k_BW:   stiffness of the BW-elements    [kN/m]
beta:   parameter of the Bouc-Wen model [-]
gamma:  parameter of the Bouc-Wen model [-]
n:      parameter of the Bouc-Wen model [-]
dt:     times step                      [-]

OUTPUTS:
xi_i: end Bouc-Wen displacements                        [m]
r_i:  restoring force vector                            [kN]
K_BW: tangent stiffness matrix of the Bouc-Wen elements [-]

MADE BY:      junsebas97
BIBLIOGRAPHY: {1} A nonlinear Bayesian filter with uncertain inputs based on
                  the generalised-alpha method - Delgado Trujillo et al
              {2} Numnerical Methods for Engineers (6th edition) - Chapra
                  SC, Canale RP
%}
N_DOFs = size(alpha, 1);    % number of DOFs
dv     = v_i - v_im1;       % velocity change [m/s2]

% define the dynamic model -- {1} Eq.80
g = @(xi, v) diff([0; v]).* ...                                          % [m/s]
             (1 - (abs(xi).^n).*(beta.*sign(xi.*diff([0; v])) + gamma));

% evaluate the end BW displacement with Runge-Kutta
k1   = g(xi_im1,               v_im1           );    % [m/s] -- {2} Eqs.25.40a
k2   = g(xi_im1 + (1/2)*k1*dt, v_im1 + (1/2)*dv);    % [m/s] -- {2} Eqs.25.40b
k3   = g(xi_im1 + (1/2)*k2*dt, v_im1 + (1/2)*dv);    % [m/s] -- {2} Eqs.25.40c
k4   = g(xi_im1 +       k3*dt, v_im1 +       dv);    % [m/s] -- {2} Eqs.25.40d

xi_i = xi_im1 + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*dt;    % [m]   -- {2} Eq.25.40

% compute the restoring forces
r_i = alpha.*k.*diff([0; u_i]) + ...    % spring forces [kN] -- {1} Eq.79
      (1 - alpha).*k.*xi_i;
r_i = r_i - [r_i(2:end); 0];            % DOFs forces   [kN]

% calculate the tangent stiffness matrix of the Bouc-Wen elements
dxi_du = 1 - (abs(xi_i).^n).*...                            % tangent of BW
             (beta.*sign(xi_i.*diff([0; v_i])) + gamma);    % displacement [-]
                                                            % from -- {1} Eq.80

kT_BW = alpha.*k + (1 - alpha).*k.*dxi_du;                  % tangent elemental
                                                            % stiffness [kN/m]

K_BW       = zeros(N_DOFs, N_DOFs);                         % tangent stiffness
K_BW(1, 1) = kT_BW(1);                                      % matrix
for i = 2:N_DOFs
    idx            = [i - 1, i];
    K_BW(idx, idx) = K_BW(idx, idx) + [ kT_BW(i), -kT_BW(i);
                                       -kT_BW(i),  kT_BW(i)];
end
end