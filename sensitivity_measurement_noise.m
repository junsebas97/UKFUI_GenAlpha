%{
This code implements the nonlinear Bayesian filter with uncertain inputs.
It compares the filter's performance with different noise levels.

MADE BY:      junsebas97
BIBLIOGRAPHY: A nonlinear Bayesian filter with uncertain inputs based on
              the generalised-alpha method - Delgado Trujillo et al
%}
close all; clear all; clc
rng(1234)

%% INPUTS:
% system:
m         = [4; 4; 4; 3; 3; 3; 3; 3];    % mass                 [1e3 kg]
c         = [1; 1; 1; 2; 2; 2; 2; 2];    % damping              [kN-s/m]
k         = 15*ones(8, 1);               % stiffness            [kN/m]
alpha_BW  = 0.5*ones(8, 1);              % stiffness ratio      [-]
beta_BW   = 60*ones(8, 1);               % BW model parameter   [-]
gamma_BW  = -40*ones(8, 1);              % BW model parameter   [-]
n_BW      = 2*ones(8, 1);                % BW model parameter   [-]
x0        = zeros(24, 1);                % initial condition    [m], [m/s], [m]
tolerance = 1e-3;                        % iterations tolerance [-]
max_iter  = 20;                          % maximum iterations   [-]

% force:
ElCentro = load('el_centro.txt');
t        = ElCentro(:, 1);           % times                             [s]
u_u      = 9.806*ElCentro(:, 2)';    % uncertain input                   [kN]
u_k      = zeros(1, size(t, 1));     % known input                       [kN]
eta_u    = -m;                       % uncertain input influence on f(t) [1e3kg]
phi_u    = zeros(8, 1);              % uncertain input influence on p(t) [-]
eta_k    = zeros(8, 1);              % known input influence on f(t)     [-]
phi_k    = eta_k;                    % known input influence on p(t)     [-]

% measurements:
noise_ratio = [0.05, 0.1, 0.25, 0.5];        % noise2signal ratio            [-]

Su = sparse([1, 2], [1, 8], [1, 1], 2, 8);   % displacement selection matrix [-]
Sv = sparse(     1,      6,      1, 1, 8);   % velocity selection matrix     [-]
Sa = sparse(     1,      3,      1, 1, 8);   % acceleration selection matrix [-]

% filter:
m_fltr     = 3*ones(8, 1);     % mass                     [1e3kg]
alpha_fltr = 0.4*ones(8, 1);   % stiffness ratio          [-]

kappa      = 1;                % spread factor UT         [-]
alpha      = 1;                % spread factor UT         [-]
beta       = 0;                % non-Gaussian factor UT   [-]

mu_0 = (1e-1);                 % initial input mean       [m/s2]
Pu_0 = (1e-1);                 % initial input covariance [m2/s4]
E    = (1e-1);                 % random walk covariance   [m2/s4]

mx_0 = (1e-1)*randn(24, 1);    % initial state mean       [m],  [m/s],   [m]
Px_0 = (1e-1)*eye(24);         % initial state covariance [m2], [m2/s2], [m2]
Q    = (1e-4)*eye(24);         % process noise covariance [m2], [m2/s2], [m2]

R    = (1e-3)*eye(4);          % measurement noise covar  [m2], [m2/s2], [m2/s4]

%% PARAMETERS:
Nt     = size(t, 1);            % number of time steps
N_DOFs = size(m, 1);            % number of DOFs
nx     = size(mx_0, 1);         % dimensionality of the state
nu     = size(mu_0, 1);         % dimensionality of the uncertain input

dt     = t(2) - t(1);           % time step [s]

% define the system parameters
M       = sparse(diag(m));           % mass matrix    [1e3 kg]
C       = sparse(N_DOFs, N_DOFs);    % damping matrix [kN-s/m]
C(1, 1) = c(1);

for i = 2:N_DOFs
    C([i - 1, i], [i - 1, i]) = C([i - 1, i], [i - 1, i]) + [ c(i), -c(i);
                                                             -c(i),  c(i)];
end

M_fltr    = sparse(diag(m_fltr));  % mass matrix for filtering         [1e3 kg]
etau_fltr = -m_fltr;               % uncertain influence for filtering [1e3 kg]

theta = {x0; []; C; k; []; beta_BW; gamma_BW; n_BW; dt; tolerance; max_iter;
         Su; Sv; Sa; u_k; eta_k; phi_k; u_u; []; phi_u};

%% MAIN:
N_analy = numel(noise_ratio);
y       = cell(N_analy, 1);
mx      = cell(N_analy, 1);
mu      = cell(N_analy, 1);
Px      = cell(N_analy, 1);
Pu      = cell(N_analy, 1);
xRRMSE  = NaN(nx, N_analy);
uRRMSE  = NaN(nu, N_analy);

for i = 1:N_analy
    % in each analysis,
    % 1) create target data (measurements and states)
    rng(1234)
    theta{2}  = M;
    theta{5}  = alpha_BW;
    theta{19} = eta_u;
    [x, y{i}] = get_data(noise_ratio(i), theta);

    % 2) apply the proposed filter
    theta{2}                     = M_fltr;
    theta{5}                     = alpha_fltr;
    theta{19}                    = etau_fltr;
    [mx{i}, Px{i}, mu{i}, Pu{i}] = this_filter(y{i}, mx_0, Px_0, mu_0, Pu_0, ...
                                               E, Q, R, kappa, alpha, beta,  ...
                                               theta);
    % 3) evaluate the RRMSE -- Eq.83
    xRRMSE(:, i) = rmse( mx{i},   x, 2)./max(  x, [], 2);    % [m], [m/s], [m]
    uRRMSE(:, i) = rmse(mu{i},  u_u, 2)./max(u_u, [], 2);    % [m/s2]

    % 4) plot the filter's estimation
    plot_estimation(  x, y{i}, t, mx{i}, Px{i}, {Su; Sv; sparse(0, N_DOFs)})
    plot_estimation(u_u,   [], t, mu{i}, Pu{i}, {        sparse(0,     nu)})
end

%{
NOTE:
The units are as follows
    x and mx ---> [m],  [m/s],   [m]
    Px       ---> [m2], [m2/s2], [m2]
    y        ---> [m],  [m/s],   [m/s2]
    mu       ---> [m/s2]
    Pu       ---> [m2/s4]
%}

%% REPORT:
compare_estimation(  x, y{end}, t, mx, {Su; Sv; sparse(0, N_DOFs)}, ["5%", "10%", "25%", "50%"])
compare_estimation(u_u,     [], t, mu, {        sparse(0,     nu)}, ["5%", "10%", "25%", "50%"])

print_RMSE(xRRMSE, ["5%", "10%", "25%", "50%"])
print_RMSE(uRRMSE, ["5%", "10%", "25%", "50%"])