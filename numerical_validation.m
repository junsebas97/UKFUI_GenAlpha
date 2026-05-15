%{
This code implements the nonlinear Bayesian filter with uncertain inputs,
and compares it agains the UKF-UI of Lei Y et al and the UKF-UL of
Delgado-Trujillo et al. It estimates the states and the uncertain inputs of the
system from noisy measurements

MADE BY:      junsebas97
BIBLIOGRAPHY: A nonlinear Bayesian filter with uncertain inputs based on
              the generalised-alpha method - Delgado Trujillo et al
%}
close all; clear all; clc
rng(1234)

%% INPUTS:
% system:
m         = [4, 4, 4, 3, 3, 3, 3, 3];    % mass                 [1e3 kg]
c         = [1, 1, 1, 2, 2, 2, 2, 2];    % damping              [kN-s/m]
k         = 15*ones(8, 1);               % stiffness            [kN/m]
alpha_BW  = 0.5*ones(8, 1);              % stiffness ratio      [-]
beta_BW   = 60*ones(8, 1);               % BW model parameter   [-]
gamma_BW  = -40*ones(8, 1);              % BW model parameter   [-]
n_BW      = 2*ones(8, 1);                % BW model parameter   [-]
x0        = zeros(24, 1);                % initial condition    [m], [m/s], [m]
tolerance = 1e-3;                        % iterations tolerance [-]
max_iter  = 20;                          % maximum iterations   [-]

% force:
t     = (0:0.05:60)';                          % times                      [s]
eta_k = sparse(     5,      1,      1, 8, 1);  % known influence matrix     [-]
eta_u = sparse([2, 7], [1, 2], [1, 1], 8, 2);  % uncertain influence matrix [-]
phi_k = eta_k;                                 % known influence matrix     [-]
phi_u = eta_u;                                 % uncertain influence matrix [-]
u_k   = 0.75*tanh(0.1*t');                     % known input                [kN]
u_u   = [randn(1, size(t, 1));                 % uncertain inputs           [kN]
         -0.04*(10 + t'.*sin(t'))];

% measurements:
noise_ratio = 0.05;                                    % noise2signal ratio  [-]

Su = sparse(   [1, 2],    [1, 6],    [1, 1], 2, 8);    % displacement select [-]
Sv = sparse(        1,         7,         1, 1, 8);    % velocity select     [-]
Sa = sparse([1, 2, 3], [2, 3, 8], [1, 1, 1], 3, 8);    % acceleration select [-]

% filter:
m_fltr     = 3*ones(8, 1);     % mass                      [1e3kg]
alpha_fltr = 0.4*ones(8, 1);   % stiffness ratio           [-]

kappa      = 1;                % spread factor UT          [-]
alpha      = 1;                % spread factor UT          [-]
beta       = 0;                % non-Gaussian factor UT    [-]

mu_0 = (1e-1)*ones(2, 1);      % initial input mean       [kN]
Pu_0 = (1e-1)*eye(2);          % initial input covariance [kN2]
E    = (1e-1)*eye(2);          % random walk covariance   [kN2]

mx_0 = (1e-1)*randn(24, 1);    % initial state mean       [m],  [m/s],   [m]
Px_0 = (1e-1)*eye(24);         % initial state covariance [m2], [m2/s2], [m2]
Q    = (1e-4)*eye(24);         % process noise covariance [m2], [m2/s2], [m2]

R    = (1e-3)*eye(6);          % measurement noise covar  [m2], [m2/s2], [m2/s4]

%% PARAMETERS:
Nt     = size(t, 1);            % number of time steps
N_DOFs = size(m, 1);            % number of DOFs
nx     = size(mx_0, 1);         % dimensionality of the state
nu     = size(mu_0, 1);         % dimensionality of the uncertain input

dt     = t(2) - t(1);           % time step [s]

% define the system parameters
M       = sparse(diag(m));           % mass matrix for data      [1e3 kg]
M_fltr  = sparse(diag(m_fltr));      % mass matrix for filtering [1e3 kg]
C       = sparse(N_DOFs, N_DOFs);    % damping matrix            [kN-s/m]
C(1, 1) = c(1);

for i = 2:N_DOFs
    C([i - 1, i], [i - 1, i]) = C([i - 1, i], [i - 1, i]) + [ c(i), -c(i);
                                                             -c(i),  c(i)];
end

theta = {x0; M; C; k; alpha_BW; beta_BW; gamma_BW; n_BW; dt; tolerance; 
         max_iter; Su; Sv; Sa; u_k; eta_k; phi_k; u_u; eta_u; phi_u};

%% MAIN:
% create target data
[x, y] = get_data(noise_ratio, theta);    % states       [m], [m/s], [m]
                                          % measurements [m], [m/s], [m/s2]

% get the estimations of the proposed filter, the UKF-UI, and the UKF-UL
theta{2} = M_fltr;
theta{5} = alpha_fltr;

[mx, Px, mu, Pu] = Ga_UKFUI(y, mx_0, Px_0, mu_0, Pu_0, E, Q, R, kappa, ...
                            alpha, beta, theta);

[mx_UKFUL, Px_UKFUL, mfu_UKFUL, Pfu_UKFUL] = UKF_UL(y, mx_0, Px_0, mu_0,  ...
                                                    Pu_0, E, Q, R, kappa, ...
                                                    alpha, beta, theta);

[mx_UKFUI, Px_UKFUI, fu_UKFUI] = UKF_UI(y, mx_0, Px_0, mu_0, Q, R, kappa, ...
                                       alpha, beta, theta);
%{
NOTE:
The units are as follows
    mx, mx_UKFUL, and mx_UKFUI  ---> [m],  [m/s],   [m]
    Px, Px_UKFUL, and Px_UKFUI  ---> [m2], [m2/s2], [m2]
    mu, mfu_UKFUL, and fu_UKFUI ---> [kN]
    Pu and Pfu_UKFUL            ---> [kN2]
%}

% evaluate the RRMSE of the estimations
xRRMSE_GaUKFUI = rmse(      mx, x, 2)./max(x, [], 2);        % [m], [m/s], [m]
xRRMSE_UKFUL   = rmse(mx_UKFUL, x, 2)./max(x, [], 2);        % -- Eq.83
xRRMSE_UKFUI   = rmse(mx_UKFUI, x, 2)./max(x, [], 2);

uRRMSE_GaUKFUI = rmse(       mu, u_u, 2)./max(u_u, [], 2);   % [kN] -- Eq.83
uRRMSE_UKFUL   = rmse(mfu_UKFUL, u_u, 2)./max(u_u, [], 2);
uRRMSE_UKFUI   = rmse( fu_UKFUI, u_u, 2)./max(u_u, [], 2);

%% REPORT:
compare_estimation(x, y, t, {mx, mx_UKFUL, mx_UKFUI}, {Su; Sv; sparse(0, N_DOFs)}, ...
                   ["G\alpha-UKFUI", "UKF-UL", "UKF-UI"])
compare_estimation(u_u, [], t, {mu, mfu_UKFUL, fu_UKFUI}, {sparse(0, nu)}, ...
                   ["G\alpha-UKFUI", "UKF-UL", "UKF-UI"])

plot_estimation(  x,  y, t, mx, Px, {Su; Sv; sparse(0, N_DOFs)})
plot_estimation(u_u, [], t, mu, Pu, {        sparse(0,     nu)})

print_RMSE([xRRMSE_GaUKFUI, xRRMSE_UKFUL, xRRMSE_UKFUI], ...
           [    "Ga-UKFUI",     "UKF-UL",     "UKF-UI"])
print_RMSE([uRRMSE_GaUKFUI, uRRMSE_UKFUL, uRRMSE_UKFUI], ...
           [    "Ga-UKFUI",     "UKF-UL",     "UKF-UI"])