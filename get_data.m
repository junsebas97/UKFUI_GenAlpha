function [x, y] = get_data(noise_ratio, theta)
%{
this function creates the data. It returns the noisy measurements, and the
true states. It uses the generalised-alpha method with rho_infty = 0.5 for
the time integration and the Bouc-Wen model for the restoring force.

INPUTS:
noise_ratio: noise-to-signal ratio of the measurements [-]
theta:       system parameters. They must organized as
   ---                                                        ---
   | initial condition                          [m], [m/s], [m] |
   | mass matrix                                [1e3 kg]        |
   | damping matrix                             [kN-s/m]        |
   | stiffnesses                                [kN/m]          |
   | post-yield stiffness ratio                 [-]             |
   | parameter of the Bouc-Wen model            [-]             |
   | parameter of the Bouc-Wen model            [-]             |
   | parameter of the Bouc-Wen model            [-]             |
   | high energy dissipation                    [-]             |
   | time step                                  [s]             |
   | tolerance Newton-Raphson                   [-]             |
   | maximum iterations Newton-Raphson          [-]             |
   | displacement selection matrix              [-]             |
   | velocity selection matrix                  [-]             |
   | acceleration selection matrix              [-]             |
   | known inputs                               [kN] or [m/s2]  |
   | known inputs influence matrix on f(t)      [-]  or [1e3kg] |
   | known inputs influence matrix on p(t)      [-]             |
   | uncertain inputs                           [kN] or [m/s2]  |
   | influence matrix uncertain inputs on f(t)  [-]  or [1e3kg] |
   | uncertain inputs influence matrix on p(t)  [-]             |
   ---                                                        ---

OUTPUTS:
x: system states [m], [m/s], [m]
y: measurements  [m], [m/s], [m/s2]

MADE BY:      junsebas97
BIBLIOGRAPHY: A nonlinear Bayesian filter with uncertain inputs based on the
              generalised-alpha method - Delgado Trujillo et al
%}

%% PARAMETRIZATION:
rho_infty = 0.5;    % high energy dissipation [-]

% extract the model parameters
x0        = theta{1};     % initial condition                  [m], [m/s], [m]
M         = theta{2};     % mass matrix                        [1e3 kg]
C         = theta{3};     % damping matrix                     [kN-s/m]
k         = theta{4};     % springs' stiffnesses               [kN/m]
alpha_BW  = theta{5};     % post-yield stiffness ratios        [-]
beta_BW   = theta{6};     % parameter of the Bouc-Wen model    [-]
gamma_BW  = theta{7};     % parameter of the Bouc-Wen model    [-]
n_BW      = theta{8};     % parameter of the Bouc-Wen model    [-]
dt        = theta{9};     % time step                          [s]
tolerance = theta{10};    % tolerance Newton-Raphson           [-]
max_iter  = theta{11};    % maximum iterations Newton-Raphson  [-]
Su        = theta{12};    % displacement selection matrix      [-]
Sv        = theta{13};    % velocity selection matrix          [-]
Sa        = theta{14};    % acceleration selection matrix      [-]
u_k       = theta{15};    % known inputs                       [kN] or [m/s2]
eta_k     = theta{16};    % known influence matrix on f(t)     [-]  or [1e3kg]
phi_k     = theta{17};    % known influence matrix on p(t)     [-]
u_u       = theta{18};    % uncertain inputs                   [kN] or [m/s2]
eta_u     = theta{19};    % uncertain influence matrix on f(t) [-]  or [1e3kg]
phi_u     = theta{20};    % uncertain influence matrix on p(t) [-]

N_DOFs = size(M, 1);      % number of DOFs
Nt     = size(u_k, 2);    % number of times

% calculate the parameters of the TIA
alpha_f = rho_infty/(rho_infty + 1);            % [-] -- Eq.12
alpha_m = (2*rho_infty - 1)/(rho_infty + 1);    % [-] -- Eq.11
gamma   = 1/2 - alpha_m + alpha_f;              % [-] -- Eq.10
beta    = (1/4)*(1 - alpha_m + alpha_f)^2;      % [-] -- Eq.9

% extract the previous state -- Page 4
u0  = x0(             1:  N_DOFs);    % displacement    [m]
v0  = x0(  (N_DOFs + 1):2*N_DOFs);    % velocity        [m/s]
xi0 = x0((2*N_DOFs + 1):     end);    % BW displacement [m]

% compute the external force and the applied loads
f = eta_k*u_k + eta_u*u_u;    % [kN] -- Eq.33
p = phi_k*u_k + phi_u*u_u;    % [kN] -- Eq.34

% assess the initial restoring force vector and the initial acceleration
r0 = alpha_BW.*k.*diff([0; u0]) + (1 - alpha_BW).*k.*xi0; % [kN] -- Eq.79
r0 = r0 - [r0(2:end);  0];                                % [kN]

a0 = M\(f(:, 1) - C*v0 - r0);                             % [m/s2] -- from Eq.78

%% TIME-INTEGRATION:
% assign the initial conditions
u  = NaN(N_DOFs, Nt);      u(:, 1)  = u0;     % [m]
v  = NaN(N_DOFs, Nt);      v(:, 1)  = v0;     % [m/s]
a  = NaN(N_DOFs, Nt);      a(:, 1)  = a0;     % [m/s2]
xi = NaN(N_DOFs, Nt);      xi(:, 1) = xi0;    % [m]
r  = NaN(N_DOFs, Nt);      r(:, 1)  = r0;     % [kN]

% simulate the deterministic system to get data
for i = 2:Nt
    % for each time step,
    % 1) sample the force with a linear approximation
    f_tau = (1 - alpha_f)*f(:, i) + alpha_f*f(:, i - 1);    % [kN]

    % 2) initalize the displacements and velocities
    u(:, i) = u(:, i - 1);    % [m]
    v(:, i) = v(:, i - 1);    % [m/s]

    % 3) compute the equivalent stiffness matrix and force vector
    K_hat = M*((1 - alpha_m)/(beta*dt^2))     + ...    % [kN/m]
            C*((1 - alpha_f)*gamma/(beta*dt));

    f_hat = f_tau                                                      - ... % [kN]
            (  M*(1 - 1/(2*beta) + alpha_m/(2*beta))                   + ...
               C*(dt*(1 - alpha_f)*(1 - gamma/(2*beta))) )*a(:, i - 1) - ...
            (  C*(1 - gamma/beta + alpha_f*gamma/beta)                 - ...
               M*((1 - alpha_m)/(beta*dt))               )*v(:, i - 1) - ...
            (- M*((1 - alpha_m)/(beta*dt^2))                           - ...
               C*((1 - alpha_f)*gamma/(beta*dt))         )*u(:, i - 1) - ...
            alpha_f*r(:, i - 1);

    % 4) find the new nodal displacements with Newton-Raphson iteration
    for iter = 1:max_iter
        % 4.1) calculate the BW response
        [xi(:, i), r(:, i), K] = BW_model(xi(:, i - 1),         ...    % [m], 
                                          v(:, i - 1), v(:, i), ...    % [kN],
                                          u(:, i), k, alpha_BW, ...    % [kN/m]
                                          beta_BW, gamma_BW, n_BW, dt);

        % 4.2) assess the tangent stiffness and compute the residual force
        K_t   = K_hat + (1 - alpha_f)*K;                            % [kN/m]
        R_for = f_hat - (K_hat*u(:, i) + (1 - alpha_f)*r(:, i));    % [kN]

        % 4.3) if tolerance is not met, update the current displacements
        %      and velocities
        if norm(R_for) > tolerance
            duN     = K_t\R_for;         % increment    [m]
            u(:, i) = u(:, i) + duN;     % displacement [m]

            v(:, i) = (gamma/(beta*dt))*(u(:, i) - u(:, i - 1)) + ... % velocity
                      (1 - gamma/beta)*v(:, i - 1)              + ... % [m/s]
                      dt*(1 - gamma/(2*beta))*a(:, i - 1);
        else
            break
        end
    end

    % 5) compute the current accelerations
    a(:, i) = (1/(beta*dt^2))*(u(:, i) - u(:, i - 1)) - ...    % [m/s2]
              (1/(beta*dt))*v(:, i - 1)               - ...
              (1/(2*beta) - 1)*a(:, i - 1);
end

%% DATA CREATION:
% define the state and get the measurements
x = [u; v; xi];                % state        [m], [m/s], [m]    -- Page 4
y = [Su*u;                     % measurements [m], [m/s], [m/s2] -- Eq.32
     Sv*v
     (Sa/M)*(p - C*v - r)];

% add noise to the measurements
noise_mean = zeros(size(y, 1), 1);                      % [m],  [m/s],   [m/s2]
noise_var  = diag(noise_ratio*var(y, 0, 2));            % [m2], [m2/s2], [m2/s4]
y          = y + mvnrnd(noise_mean, noise_var, Nt)';    % [m],  [m/s],   [m/s2]
end