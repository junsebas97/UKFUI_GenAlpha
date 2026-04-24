function [mx, Px, mfu, Pfu] = UKF_UL(y, mx_0, Px_0, mfu_0, Pfu_0, E, Q, R, ...
                                     kappa, alpha, beta, theta)
%{
this function is the nonlinear Bayesian filter with uncertain loads. It
infers the states and the uncertain forces of the given system using ramdon
walks, the unscented transform, and SSTIAs.

INPUTS:
y:     measurements                              [m],  [m/s],   [m/s2]
mx_0:  initial state mean                        [m],  [m/s],   [m]
Px_0:  initial state covariance                  [m2], [m2/s2], [m2]
mfu_0: mean of the initial uncertain force       [kN]
Pfu_0: covariance of the initial uncertain force [kN2]
E:     covariance of uncertain force random walk [kN2]
Q:     covariance process noise                  [m2], [m2/s2], [m2]
R:     covariance measurement noise              [m2], [m2/s2], [m2/s4]
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
mx:  state mean                     [m],  [m/s],   [m]
Px:  state covariance               [m2], [m2/s2], [m2]
mfu: mean of uncertain forces       [kN]
Pfu: covariance of uncertain forces [kN]

MADE BY:      junsebas97
BIBLIOGRAPHY: {1} A nonlinear Bayesian filter for structural systems with
                  uncertain forces - Delgado Trujillo et al
              {2} A Time Integration Algorithm for Structural Dynamics With
                  Improved Numerical Dissipation: The Generalized-alpha
                  Method - Chung J, Hulbert G
%}
%% PARAMETRIZATION:
rho_infty = 0.8;         % high energy dissipation of the TIA    [-]

% extract and reorganise the system parameters
M         = theta{2};     % mass matrix                         [1e3 kg]
C         = theta{3};     % damping matrix                      [kN-s/m]
k         = theta{4};     % springs' stiffnesses                [kN/m]
alpha_BW  = theta{5};     % post-yield stiffness ratios         [-]
beta_BW   = theta{6};     % parameter of the Bouc-Wen model     [-]
gamma_BW  = theta{7};     % parameter of the Bouc-Wen model     [-]
n_BW      = theta{8};     % parameter of the Bouc-Wen model     [-]
dt        = theta{9};     % time step                           [s]
tolerance = theta{10};    % tolerance Newton-Raphson            [-]
max_iter  = theta{11};    % maximum iterations Newton-Raphson   [-]
Su        = theta{12};    % displacement selection matrix       [-]
Sv        = theta{13};    % velocity selection matrix           [-]
Sa        = theta{14};    % acceleration selection matrix       [-]
f_k       = theta{15};    % known forces                        [kN]
S_fk      = theta{16};    % collocation matrix known forces     [-]
S_fu      = theta{19};    % collocation matrix uncertain inputs [-]

theta = {M; C; k; alpha_BW; beta_BW; gamma_BW; n_BW; rho_infty; dt;
         tolerance; max_iter};

% add the acceleration components to the state variables
N_DOFs = size(M, 1);                       % number of DOFs

mx_0   = [mx_0(           (1:N_DOFs));     % displacements   [m]
          mx_0(  N_DOFs + (1:N_DOFs));     % velocities      [m/s]
          zeros(N_DOFs, 1);                % accelerations   [m/s2]
          mx_0(2*N_DOFs + (1:N_DOFs))];    % BW displacement [m]

P_aux                                                   = zeros(4*N_DOFs);
P_aux(           1:(2*N_DOFs),            1:(2*N_DOFs)) = Px_0(      1:(2*N_DOFs),       1:(2*N_DOFs));
P_aux(           1:(2*N_DOFs),      (3*N_DOFs + 1):end) = Px_0(      1:(2*N_DOFs), (2*N_DOFs + 1):end);
P_aux(     (3*N_DOFs + 1):end,            1:(2*N_DOFs)) = Px_0((2*N_DOFs + 1):end,       1:(2*N_DOFs));
P_aux(     (3*N_DOFs + 1):end,      (3*N_DOFs + 1):end) = Px_0((2*N_DOFs + 1):end, (2*N_DOFs + 1):end);
P_aux((2*N_DOFs + 1):3*N_DOFs, (2*N_DOFs + 1):3*N_DOFs) = mean(Px_0(Px_0 > 0), 'all')*eye(N_DOFs, N_DOFs);
Px_0                                                    = P_aux;

Q_aux                                                   = zeros(4*N_DOFs);
Q_aux(           1:(2*N_DOFs),            1:(2*N_DOFs)) = Q(      1:(2*N_DOFs),       1:(2*N_DOFs));
Q_aux(           1:(2*N_DOFs),      (3*N_DOFs + 1):end) = Q(      1:(2*N_DOFs), (2*N_DOFs + 1):end);
Q_aux(     (3*N_DOFs + 1):end,            1:(2*N_DOFs)) = Q((2*N_DOFs + 1):end,       1:(2*N_DOFs));
Q_aux(     (3*N_DOFs + 1):end,      (3*N_DOFs + 1):end) = Q((2*N_DOFs + 1):end, (2*N_DOFs + 1):end);
Q_aux((2*N_DOFs + 1):3*N_DOFs, (2*N_DOFs + 1):3*N_DOFs) = mean(Q(Q > 0), 'all')*eye(N_DOFs, N_DOFs);
Q                                                       = Q_aux;

% get lengths
Nt     = size(f_k, 2);          % number of time steps
N_Su   = size(Su, 1);           % number displacement measurements
N_Sv   = size(Sv, 1);           % number velocities measurements
N_Sa   = size(Sa, 1);           % number acceleration measurements
nx     = size(mx_0, 1);         % dimensionality of the state
nfu    = size(mfu_0, 1);        % dimensionality of the uncerain force
nz     = nx + nfu;              % augmented dimensionality -- {1} Page 4
Ny     = N_Su + N_Sv + N_Sa;    % total number of measurements

% evaluate the known force at the sampling time
alpha_f       = rho_infty/(rho_infty + 1);          % force factor [-] -- {2} Eq.25
f_k(:, 2:end) = alpha_f*f_k(:, (2:end) - 1) + ...   % known force  [kN]
                (1 - alpha_f)*f_k(:, (2:end));

% construct the measurement matrix
H                                                = zeros(Ny, nx);    % -- Eq.12
H(              (1:N_Su),            (1:N_DOFs)) = Su;               % [-]
H(       N_Su + (1:N_Sv),   N_DOFs + (1:N_DOFs)) = Sv;               % [-]
H(N_Su + N_Sv + (1:N_Sa), 2*N_DOFs + (1:N_DOFs)) = Sa;               % [-]

% compute the weights of the UT
lambda    = alpha^2*(nz + kappa) - nz;                       % [-] -- {1} Page 3
Wm        = NaN(2*nz + 1, 1);
Wc        = NaN(1, 2*nz + 1);
Wm(1)     = lambda/(nz + lambda);                            % [-] -- {1} Eq.7
Wc(1)     = (lambda/(nz + lambda)) + (1 - alpha^2 + beta);   % [-] -- {1} Eq.8
Wm(2:end) = 1/(2*(nz + lambda));                             % [-] -- {1} Eq.9
Wc(2:end) = 1/(2*(nz + lambda));                             % [-] -- {1} Eq.9

%% FILTERING:
% apply the filter to estimate the sytem state and the uncertain force
mx  = NaN( nx,      Nt);       mx( :,    1) = mx_0;
mfu = NaN(nfu,      Nt);       mfu(:,    1) = mfu_0;
Px  = NaN( nx,  nx, Nt);       Px( :, :, 1) = Px_0;
Pfu = NaN(nfu, nfu, Nt);       Pfu(:, :, 1) = Pfu_0;

for i = 2:Nt
    % for each time step:
    % 1) predict the input
    mfu_iim1 = mfu(:,    i - 1);        % mean       [kN]  -- {1} Eq.34
    Pfu_iim1 = Pfu(:, :, i - 1) + E;    % covariance [kN2] -- {1} Eq.35

    % 2) predict the state and compute its covariances
    % 2.1) create the sigma points of the UT
    mz      = [mx(:, i - 1); mfu_iim1];                         %      {1} Eq.20
    sqrt_Pz = sqrt(nz + lambda)*...                             % from {1} Eq.21
              chol(blkdiag(Px(:, :, i - 1), Pfu_iim1), 'lower');

    Z                     = NaN(nz, 2*nz + 1);
    Z(:,               1) = mz;                   % {1} Eq.1
    Z(:, 1 +      (1:nz)) = mz + sqrt_Pz;         % {1} Eq.2
    Z(:, 1 + nz + (1:nz)) = mz - sqrt_Pz;         % {1} Eq.3

    X = Z(      (1:nx), :);         % [m], [m/s], [m/s2], [m] -- {1} Page 4
    U = Z(nx + (1:nfu), :);         % [kN]                    -- {1} Page 4
    F = S_fk*f_k(:, i) + S_fu*U;    % [kN]                    -- {1} Page 4

    % 2.2) evaluate the current state by propagating the sigma points
    g_XF = NaN(nx, 2*nz + 1);
    for j = 1:(2*nz + 1)
        g_XF(:, j) = g(X(:, j), F(:, j), theta);   % [m], [m/s], [m/s2], [m]
    end

    % 2.3) compute the statistics of the current state
    mx_iim1 = g_XF*Wm;                        % mean       [m],     -- {1} Eq.36
                                              %            [m/s], 
                                              %            [m/s2],
                                              %            [m]
    g_dist  = g_XF - mx_iim1;
    Px_iim1 = (Wc.*g_dist)*g_dist' + Q;       % covariance [m2],    -- {1} Eq.37
                                              %            [m2/s2],
                                              %            [m2/s4],
                                              %            [m2]
    
    % 3) update the state and uncertain force -- {1} Eqs.44 to 52
    my           = H*mx_iim1;                     % [m],    [m/s],    [m/s2]
    Py           = H*Px_iim1*H' + R;              % [m2],   [m2/s2],  [m2/s4]
    P_fux        = (Wc.*(U - mfu_iim1))*g_dist';  % [kN-m], [kN-m/s], [kN-m/s2], [kN-m]

    Kf           = (P_fux*H')/Py;
    mfu(:,    i) = mfu_iim1 + Kf*(y(:, i) - my);  % [kN]
    Pf_aux       = Pfu_iim1 - Kf*Py*Kf';          % [kN2]
    Pfu(:, :, i) = (1/2)*(Pf_aux + Pf_aux');

    Kx           = (Px_iim1*H')/Py;
    mx(:,    i)  = mx_iim1 + Kx*(y(:, i)  - my);  % [m],  [m/s],   [m/s2],  [m]
    Px_aux       = Px_iim1 - Kx*Py*Kx';           % [m2], [m2/s2], [m2/s4], [m2]
    Px(:, :, i)  = (1/2)*(Px_aux + Px_aux');
end

%% FORMATING:
% remove the acceleration from the estimations
mx(2*N_DOFs + (1:N_DOFs),                        :) = [];   % mean       [m/s2]
Px(2*N_DOFs + (1:N_DOFs),                     :, :) = [];   % covariance [m2/s4]
Px(                    :, 2*N_DOFs + (1:N_DOFs), :) = [];
end

function [x_i] = g(x_im1, f_i, theta)
%{
this function is the SSTIA-built process model. It computes the state of the
given system with the previous state and the current force using the
generalised-alpha method. The system uses the Bouc-Wen model for the
restoring force.

INPUTS:
x_im1: previous system state  [m], [m/s], [m/s2], [m]
f_i:   current force          [kN]
theta: system parameters.

OUTPUTs:
x_i: current system state  [m], [m/s], [m/s2], [m]

MADE BY:      junsebas97
BIBLIOGRAPHY: {1} A nonlinear Bayesian filter for structural systems with
                  uncertain forces - Delgado Trujillo et al
              {2} A Time Integration Algorithm for Structural Dynamics With
                  Improved Numerical Dissipation: The Generalized-alpha
                  Method - Chung J, Hulbert G
%}

% get the system properties
M         = theta{1};      % mass matrix                       [1e3 kg]
C         = theta{2};      % damping matrix                    [kN-s/m]
k         = theta{3};      % stiffnesses                       [kN/m]
alpha_BW  = theta{4};      % post-yield stiffness ratio        [-]
beta_BW   = theta{5};      % parameter of the Bouc-Wen model   [-]
gamma_BW  = theta{6};      % parameter of the Bouc-Wen model   [-]
n_BW      = theta{7};      % parameter of the Bouc-Wen model   [-]
rho_infty = theta{8};      % high energy dissipation           [-]
dt        = theta{9};      % time step                         [s]
tolerance = theta{10};     % tolerance Newton-Raphson          [-]
max_iter  = theta{11};     % maximum iterations Newton-Raphson [-]

N_DOFs    = size(M, 1);    % numbers of DOFs

% calculate the parameters of the generalised-alpha method
alpha_f = rho_infty/(rho_infty + 1);            % [-] -- {2} Eq.25
alpha_m = (2*rho_infty - 1)/(rho_infty + 1);    % [-] -- {2} Eq.25
gamma   = 1/2 - alpha_m + alpha_f;              % [-] -- {2} Eq.17
beta    = (1/4)*(1 - alpha_m + alpha_f)^2;      % [-] -- {2} Eq.20

% extract the previous state -- {1} Page 3
u_tim1  = x_im1(           (1:N_DOFs));    % displacement      [m]
v_tim1  = x_im1(  N_DOFs + (1:N_DOFs));    % velocity          [m/s]
a_tim1  = x_im1(2*N_DOFs + (1:N_DOFs));    % acceleration      [m/s2]
xi_tim1 = x_im1(3*N_DOFs + (1:N_DOFs));    % BW displacements  [m]

% assess the previous restoring force
r_tim1 = alpha_BW.*k.*diff([0; u_tim1]) + ...    % springs force [kN] 
         (1 - alpha_BW).*k.*xi_tim1;             % -- {1} Eq.54
r_tim1 = r_tim1 - [r_tim1(2:end); 0];            % DOFs force    [kN]

% initalize the displacements and velocities
u_ti = u_tim1;    % [m]
v_ti = v_tim1;    % [m/s]

% compute the equivalent stiffness matrix and force vector
K_hat = M*((1 - alpha_m)/(beta*dt^2))     +                      ... % [kN/m]
        C*((1 - alpha_f)*gamma/(beta*dt));

f_hat = f_i                                                    - ... % [kN]
        (  M*(1 - 1/(2*beta) + alpha_m/(2*beta))               + ...
           C*(dt*(1 - alpha_f)*(1 - gamma/(2*beta)))  )*a_tim1 - ...
        (  C*(1 - gamma/beta + alpha_f*gamma/beta)             - ...
           M*((1 - alpha_m)/(beta*dt))                )*v_tim1 - ...
        (- M*((1 - alpha_m)/(beta*dt^2))                       - ...
           C*((1 - alpha_f)*gamma/(beta*dt))          )*u_tim1 - ...
        alpha_f*r_tim1;

% find the new nodal displacements with Newton-Raphson iteration
for iter = 1:max_iter
    % in each iteration,
    % 1) calculate the BW response
    [xi_ti, r_ti, dr_du] = BW_model(xi_tim1, v_tim1, v_ti, u_ti, ...   % [m],
                                    k, alpha_BW, beta_BW,        ...   % [kN],
                                    gamma_BW, n_BW, dt);               % [kN/m],

    % 2) assess the tangent stiffness and compute the residual force
    K_t   = K_hat + (1 - alpha_f)*dr_du;                  % [kN/m]
    R_for = f_hat - (K_hat*u_ti + (1 - alpha_f)*r_ti);    % [kN]

    % 3) if tolerance is not met, update the current displacements and
    %    velocities
    if norm(R_for) > tolerance
        duN  = K_t\R_for;                                 % increment    [m]
        u_ti = u_ti + duN;                                % displacement [m]
        vt_i = (gamma/(beta*dt))*(u_ti - u_tim1) + ...    % velocity     [m/s]
               (1 - gamma/beta)*v_tim1           + ... 
               dt*(1 - gamma/(2*beta))*a_tim1;
    else
        break
    end
end

% compute the current accelerations and form the current state vector
at_i = (1/(beta*dt^2))*(u_ti - u_tim1) - ...    % acceleration [m/s2] 
       (1/(beta*dt))*v_tim1            - ...
       (1/(2*beta) - 1)*a_tim1;

x_i  = [u_ti; vt_i; at_i; xi_ti];               % state [m], [m/s], [m/s2], [m]
end                                             %  -- {1} Page 3