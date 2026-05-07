%{
This code implements the nonlinear Bayesian filter with uncertain inputs.
It compares the filter's performance with different noise levels.

MADE BY:      junsebas97
BIBLIOGRAPHY: A nonlinear Bayesian filter with uncertain inputs based on
              the generalised-alpha method - Delgado Trujillo et al

% NOTE:
The measurement configurations were chosen as follows:

yvar = ["u1", "u2", "u3", "u4", "u5", "u6", "u7", "u8", "none";
        "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "none";
        "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "none"];
rng(0)
for i = 1:5
    yvar(randperm(27, 4))
end
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
noise_ratio = [0.05, 0.1, 0.25, 0.5];          % noise2signal ratio [-]

Su = {sparse(     1,      3,      1, 1, 8);    % displacement selection
      sparse([1, 2], [2, 8], [1, 1], 2, 8);    % matrixes [-]
      sparse([1, 2], [2, 5], [1, 1], 2, 8);
      sparse(0, 8);
      sparse(     1,      2,      1, 1, 8)};
Sv = {sparse(     1,      5,      1, 1, 8);    % velocity selection
      sparse(0, 8);                            % matrixes [-]
      sparse(0, 8);
      sparse(     1,      7,      1, 1, 8);
      sparse(0, 8)};
Sa = {sparse([1, 2], [1, 6], [1, 1], 2, 8);    % acceleration selection
      sparse(     1,      8,      1, 1, 8);    % matrixes [-]
      sparse(     1,      7,      1, 1, 8);
      sparse([1, 2], [4, 8], [1, 1], 2, 8);
      sparse(     1,      8,      1, 1, 8)};

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

R    = {(1e-3)*eye(4);         % measurement noise covar  [m2], [m2/s2], [m2/s4]
        (1e-3)*eye(3);         %                          [m2],          [m2/s4]
        (1e-3)*eye(3);         %                          [m2],          [m2/s4]
        (1e-3)*eye(3);         %                                [m2/s2], [m2/s4]
        (1e-3)*eye(2)};        %                          [m2],          [m2/s4]

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

theta = {x0; []; C; k; []; beta_BW; gamma_BW; n_BW; dt; tolerance;
         max_iter; []; []; []; u_k; eta_k; phi_k; u_u; []; phi_u};

%% MAIN:
N_noise = numel(noise_ratio);
N_measu = numel(Su);
y       = cell(   N_measu, N_noise);
mx      = cell(   N_measu, N_noise);
mu      = cell(   N_measu, N_noise);
Px      = cell(   N_measu, N_noise);
Pu      = cell(   N_measu, N_noise);
xRRMSE  = NaN(nx, N_measu, N_noise);
uRRMSE  = NaN(nu, N_measu, N_noise);

for i = 1:N_measu
    % for each measurement case,
    % 1) incorporate the corresponding selection matrixes
    theta{12} = Su{i};
    theta{13} = Sv{i};
    theta{14} = Sa{i};

    % 2) create date with different noise levels
    theta{2}  = M;
    theta{5}  = alpha_BW;
    theta{19} = eta_u;

    for j = 1:N_noise
        rng(1234)
        [x, y{i, j}] = get_data(noise_ratio(j), theta);
    end

    % 3) incorporate the erroneous parameters and analyse all the datasets
    theta{2}  = M_fltr;
    theta{5}  = alpha_fltr;
    theta{19} = etau_fltr;

    for j = 1:N_noise
        % 3.1) apply the proposed filter
        [mx{i, j}, Px{i, j}, mu{i, j}, Pu{i, j}] = this_filter(y{i, j}, mx_0, Px_0, ...
                                                               mu_0, Pu_0, E, Q,    ...
                                                               R{i}, kappa, alpha,  ...
                                                               beta, theta);
        % 3.2) evaluate the RRMSE -- Eq.83
        uRRMSE(:, i, j) = rmse(mu{i, j}, u_u, 2)./max(u_u, [], 2);    % [m/s2]
        xRRMSE(:, i, j) = rmse(mx{i, j},   x, 2)./max(  x, [], 2);    % [m],
                                                                      % [m/s],
                                                                      % [m]
    end

    % 4) compare the filter estimations
    compare_estimation(x, y{i, end}, t, mx(i, :),         ...
                       {Su{i}; Sv{i}; sparse(0, N_DOFs)}, ...
                       ['Case ', num2str(i), '-'] + ["5%", "10%", "25%", "50%"])
    compare_estimation(u_u, [], t, mu(i, :), {sparse(0, nu)}, ...
                       ['Case ', num2str(i), '-'] + ["5%", "10%", "25%", "50%"])
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
figure
hold on
for i = 1:N_measu
    for j = 1:N_noise
        subplot(N_measu, N_noise, j + N_noise*(i - 1))
        hold on

        for k = 1:N_DOFs
            RRMSE_DOF = xRRMSE(k:N_DOFs:end, i, j);
            plot([k, k], [min(RRMSE_DOF), max(RRMSE_DOF)], 'k-')
        end

        plot(xRRMSE(           (1:N_DOFs), i, j), 'b-', 'DisplayName',   'u(t)')
        plot(xRRMSE(  N_DOFs + (1:N_DOFs), i, j), 'r-', 'DisplayName',   'v(t)')
        plot(xRRMSE(2*N_DOFs + (1:N_DOFs), i, j), 'g-', 'DisplayName', '\xi(t)')
        xlabel 'DOF'
        ylabel 'RRMSE'
        xlim([0.5,           N_DOFs + 0.5])
        ylim([  0, max(xRRMSE, [], 'all')])
        grid on
    end
end