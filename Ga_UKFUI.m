function [mx, Px, mu, Pu] = Ga_UKFUI(y, mx_0, Px_0, mu_0, Pu_0, E, Q, R, ...
                                     kappa, alpha, beta, theta)
%{
This code implements the nonlinear Bayesian filter with uncertain inputs
based on the generalised-alpha method. It estimates the states and the
uncertain inputs of the system using the unscented transform, random walks,
and the generalised-alpha method with rho_infty = 1.

INPUTS:
y:     measurements                              [m],     [m/s],   [m/s2]
mx_0:  initial state mean                        [m],     [m/s],   [m]
Px_0:  initial state covariance                  [m2],    [m2/s2], [m2]
mu_0:  mean of the initial uncertain input       [kN]  or [m/s2]
Pu_0:  covariance of the initial uncertain input [kN2] or [m2/s4]
E:     covariance of uncertain input random walk [kN2] or [m2/s4]
Q:     covariance process noise                  [m2],    [m2/s2], [m2]
R:     covariance measurement noise              [m2],    [m2/s2], [m2/s4]
kappa: spread parameter UT                       [-]
alpha: spread parameter UT                       [-]
beta:  non-Gaussianity parameter UT              [-]
theta: system parameters. They must organized as
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

OUTPUTs:
mx:  state mean                     [m],     [m/s],   [m]
Px:  state covariance               [m2],    [m2/s2], [m2]
mfu: mean of uncertain inputs       [kN]  or [m/s2]
Pfu: covariance of uncertain inputs [kN2] or [m2/s4]

MADE BY:      junsebas97
BIBLIOGRAPHY: A nonlinear Bayesian filter with uncertain inputs based on the
              generalised-alpha method - Delgado Trujillo et al
%}

nx = size(mx_0, 1);    % dimensionality of the state
nu = size(mu_0, 1);    % dimensionality of the uncertain inputs
nz = nx + nu;          % dimensionality of the augmented vector -- Page 5
ny = size(y, 1);       % dimensionality of the measurements
nt = size(y, 2);       % number of time steps

% extract the and reorganize the system parameters
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
eta_u     = theta{19};    % uncertain influence matrix on f(t) [-]  or [1e3kg]
phi_u     = theta{20};    % uncertain influence matrix on p(t) [-]

theta_h   = {M; C; k; alpha_BW; Su; Sv; Sa};             % for measurement model
theta_g   = {M; C; k; alpha_BW; beta_BW; gamma_BW; ...   % for process model
             n_BW; dt; tolerance; max_iter};

% compute the weights of the UT [-]
lambda    = alpha^2*(nz + kappa) - nz;                        % [-] -- Page 3
Wm        = NaN(2*nz + 1, 1);
Wc        = NaN(1, 2*nz + 1);
Wm(1)     = lambda/(nz + lambda);                             % [-] -- Eq.27
Wc(1)     = (lambda/(nz + lambda)) + (1 - alpha^2 + beta);    % [-] -- Eq.28
Wm(2:end) = 1/(2*(nz + lambda));                              % [-] -- Eq.29
Wc(2:end) = 1/(2*(nz + lambda));                              % [-] -- Eq.29

% estimate the sytem state and the input with the nonlinear output-only filter
mx = NaN(nx,     nt);       mx(:,    1) = mx_0;
mu = NaN(nu,     nt);       mu(:,    1) = mu_0;
Px = NaN(nx, nx, nt);       Px(:, :, 1) = Px_0;
Pu = NaN(nu, nu, nt);       Pu(:, :, 1) = Pu_0;

for i = 2:nt
    % for each time step:
    % 1) predict the input
    mu_iim1 = mu(:,    i - 1);        % mean       [kN]  or [m/s2]  -- Eq.53
    Pu_iim1 = Pu(:, :, i - 1) + E;    % covariance [kN2] or [m2/s4] -- Eq.54

    % 2) predict the state and compute the cross-covariance with the input
    % 2.1) create the sigma points
    mz      = [mx(:, i - 1); mu_iim1];                           %      -- Eq.43
    sqrt_Pz = sqrt(nz + lambda)*...                              % from -- Eq.44
              chol(blkdiag(Px(:, :, i - 1), Pu_iim1), 'lower');

    Z                     = NaN(nz, 2*nz + 1);
    Z(:,               1) = mz;                   % -- Eq.22
    Z(:, 1 +      (1:nz)) = mz + sqrt_Pz;         % -- Eq.23
    Z(:, 1 + nz + (1:nz)) = mz - sqrt_Pz;         % -- Eq.24

    X = Z(     (1:nx), :);            % [m],    [m/s], [m] -- Page 5
    U = Z(nx + (1:nu), :);            % [kN] or [m/s2]     -- Page 5
    F = eta_k*u_k(:, i) + eta_u*U;    % [kN]               -- Page 6

    % 2.2) evaluate the current state
    g_XF = NaN(nx, 2*nz + 1);
    for j = 1:(2*nz + 1)
        g_XF(:, j) = g(X(:, j), F(:, j), theta_g);    % [m], [m/s], [m]
    end

    % 2.3) perform the prediction
    mx_iim1 = g_XF*Wm;                        % mean        [m],
                                              % -- Eq.53    [m/s],
                                              %             [m]
    g_dist  = g_XF - mx_iim1;
    Px_iim1 = (Wc.*g_dist)*g_dist' + Q;       % covariance  [m2],
                                              % -- Eq.54    [m2/s2],
                                              %             [m2]
    P_ux    = (Wc.*(U - mu_iim1))*g_dist';    % cross-covar [kN-m]   or [m2/s2],
                                              % with input  [kN-m/s] or [m2/s3],
                                              % -- Eq.50    [kN-m]   or [m2/s2]

    % 3) update the predictions
    % 3.1) create the sigma points
    mz      = [mx_iim1; mu_iim1];                                %      -- Eq.57
    sqrt_Pz = sqrt(nz + lambda)* ...                             % from -- Eq.58
              chol([Px_iim1, P_ux'; P_ux, Pu_iim1], 'lower');

    Z                     = NaN(nz, 2*nz + 1);
    Z(:,               1) = mz;                   % -- Eq.22
    Z(:, 1 +      (1:nz)) = mz + sqrt_Pz;         % -- Eq.23
    Z(:, 1 + nz + (1:nz)) = mz - sqrt_Pz;         % -- Eq.24

    X = Z(     (1:nx), :);            % [m],    [m/s], [m] -- Page 7
    U = Z(nx + (1:nu), :);            % [kN] or [m/s2]     -- Page 7
    P = phi_k*u_k(:, i) + phi_u*U;    % [kN]               -- Page 7

    % 3.2) assess the measurements
    h_XP = NaN(ny, 2*nz + 1);
    for j = 1:(2*nz + 1)
        h_XP(:, j) = h(X(:, j), P(:, j), theta_h);    % [m], [m/s], [m/s2]
    end
    
    % 3.3) compute the statistics
    my     = h_XP*Wm;                        % mean        [m],
                                             % -- Eq.74    [m/s],
                                             %             [m/s2]
    h_dist = h_XP - my;
    Py     = (Wc.*h_dist)*h_dist' + R;       % covariance  [m2],
                                             % -- Eq.75    [m2/s2],
                                             %             [m2/s4]
    P_xy   = (Wc.*(X - mx_iim1))*h_dist';    % cross-covar [m2],
                                             % with state  [m2/s],
                                             % -- Eq.76    [m2/s2],
                                             %             [m2/s3],
    P_uy   = (Wc.*(U - mu_iim1))*h_dist';    % cross-covar [kN-m]    or [m2/s2],
                                             % with input  [kN-m/s]  or [m2/s3],
                                             % -- Eq.77    [kN-m/s2] or [m2/s4]

    % 3.4) perform the update
    Ku          = P_uy/Py;                        %                     -- Eq.68
    mu(:,    i) = mu_iim1 + Ku*(y(:, i) - my);    % [kN]  or [m/s2]     -- Eq.69
    Pu(:, :, i) = Pu_iim1 - Ku*Py*Ku';            % [kN2] or [m/s2]     -- Eq.70

    Kx          = P_xy/Py;                        %                     -- Eq.71
    mx(:,    i) = mx_iim1 + Kx*(y(:, i)  - my);   % [m],  [m/s],   [m]  -- Eq.72
    Px(:, :, i) = Px_iim1 - Kx*Py*Kx';            % [m2], [m2/s2], [m2] -- Eq.73
end
end

function [x_i] = g(x_im1, f_i, theta)
%{
this function is the process model constructed with the generalised-alpha
method. It computes the state of the given system with the previous state
and the current external force. The system uses the Bouc-Wen model for the
restoring force.

INPUTS:
x_im1: previous state         [m], [m/s], [m]
f_i:   current external force [kN]
theta: system parameters:
    ---                                      ---
   | mass matrix                       [1e3 kg] |
   | damping matrix                    [kN-s/m] |
   | springs' stiffnesses              [kN/m]   |
   | post-yield stiffness ratios       [-]      |
   | parameter of the Bouc-Wen model   [-]      |
   | parameter of the Bouc-Wen model   [-]      |
   | parameter of the Bouc-Wen model   [-]      |
   | time step                         [s]      |
   | tolerance Newton-Raphson          [-]      |
   | maximum iterations Newton-Raphson [-]      |
   ---                                       ---

OUTPUTs:
xt_i: current state [m], [m/s], [m]

MADE BY:      junsebas97
BIBLIOGRAPHY: A nonlinear Bayesian filter with uncertain inputs based on the
              generalised-alpha method - Delgado Trujillo et al
%}
% get the system properties
M         = theta{1};     % mass matrix                       [1e3 kg]
C         = theta{2};     % damping matrix                    [kN-s/m]
k         = theta{3};     % springs' stiffnesses              [kN/m]
alpha     = theta{4};     % post-yield stiffness ratios       [-]
beta      = theta{5};     % parameter of the Bouc-Wen model   [-]
gamma     = theta{6};     % parameter of the Bouc-Wen model   [-]
n         = theta{7};     % parameter of the Bouc-Wen model   [-]
dt        = theta{8};     % time step                         [s]
tolerance = theta{9};     % tolerance Newton-Raphson          [-]
max_iter  = theta{10};    % maximum iterations Newton-Raphson [-]

N_DOFs    = size(M, 1);    % numbers of DOFs

% extract the previous state -- Page 4
u_tim1  = x_im1(           (1:N_DOFs));    % displacement            [m]
v_tim1  = x_im1(  N_DOFs + (1:N_DOFs));    % velocity                [m/s]
xi_tim1 = x_im1(2*N_DOFs + (1:N_DOFs));    % hysteretic displacement [m]

% assess the previous restoring force
r_tim1 = alpha.*k.*diff([0; u_tim1]) + ...    % springs force [kN] -- Eq.79
         (1 - alpha).*k.*xi_tim1;
r_tim1 = r_tim1 - [r_tim1(2:end); 0];         % DOFs force    [kN]

% initalize the displacements and velocities
u_ti = u_tim1;    % [m]
v_ti = v_tim1;    % [m/s]

% compute the equivalent stiffness matrix and force vector
K_hat = (2/(dt^2))*M + (1/dt)*C;                  % [kN/m] -- Eq.20
f_hat = f_i + (2/dt)*M*v_tim1            + ...    % [kN]   -- Eq.21
        ((2/(dt^2))*M + (1/dt)*C)*u_tim1 - ...
        (1/2)*r_tim1;

% find the new nodal displacements with Newton-Raphson iteration
for iter = 1:max_iter
    % 1) calculate the BW response
    [xi_ti, r_ti, dr_du] = BW_model(xi_tim1, v_tim1, v_ti, u_ti, ...    % [m],
                                    k, alpha, beta, gamma, n, dt);      % [kN],
                                                                        % [kN/m]
    % 2) assess the tangent stiffness and compute the residual force
    K_t   = K_hat + (1/2)*dr_du;                  % [kN/m] -- Eq.20
    R_for = f_hat - (K_hat*u_ti + (1/2)*r_ti);    % [kN]

    % 3) if tolerance is not met, update the current displacements and
    %    velocities
    if norm(R_for) > tolerance
        duN  = K_t\R_for;                          % increment    [m]
        u_ti = u_ti + duN;                         % displacement [m]
        v_i  = (2/dt)*(u_ti - u_tim1) - v_tim1;    % velocity     [m/s] -- Eq.21
    else
        break
    end
end

% form the current state vector
x_i = [u_ti; v_i; xi_ti];    % [m], [m/s], [m] -- Page 4
end

function [y_i] = h(x_i, p_i, theta)
%{
this function is the measurement model of the proposed filter. It computes
the measurements of the given system with the current state and external
force. The system uses the Bouc-Wen model for the restoring force.

INPUTS:
x_i:   state         [m], [m/s], [m]
p_i:   applied loads [kN]
theta: parameters of the measurement model:
    ---                                   ---
   | mass matrix                   [1e3 kg] |
   | damping matrix                [kN-s/m] |
   | springs' stiffnesses          [kN/m]   |
   | post-yield stiffness ratios   [-]      |
   | displacement selection matrix [-]      |
   | velocity selection matrix     [-]      |
   | acceleration selection matrix [-]      |
   ---                                    ---

OUTPUTs:
y_i: measurements [m], [m/s], [m/s2]

MADE BY:      junsebas97
BIBLIOGRAPHY: A nonlinear Bayesian filter with uncertain inputs based on the
              generalised-alpha method - Delgado Trujillo et al
%}

% get the model parameters
M      = theta{1};    % mass matrix                   [1e3 kg]
C      = theta{2};    % damping matrix                [kN-s/m]
k      = theta{3};    % spring stiffness              [kN/m]
alpha  = theta{4};    % post-yield stiffness ratio    [-]
Su     = theta{5};    % displacement selection matrix [-]
Sv     = theta{6};    % velocity selection matrix     [-]
Sa     = theta{7};    % acceleration selection matrix [-]

N_DOFs = size(M, 1);    % numbers of DOFs

% extract the state -- Page 4
u_ti  = x_i(           (1:N_DOFs));    % displacement            [m]
v_ti  = x_i(  N_DOFs + (1:N_DOFs));    % velocity                [m/s]
xi_ti = x_i(2*N_DOFs + (1:N_DOFs));    % hysteretic displacement [m]

% compute the restoring force
r_ti = alpha.*k.*diff([0; u_ti]) + ...   % springs force [kN] -- Eq.79
       (1 - alpha).*k.*xi_ti;
r_ti = r_ti - [r_ti(2:end); 0];          % DOFs force    [kN]

% calculate the measurements
y_i = [Su*u_ti;                          % [m],   -- Eq.32
       Sv*v_ti;                          % [m/s],
       (Sa/M)*(p_i - C*v_ti - r_ti)];    % [m/s2]
end