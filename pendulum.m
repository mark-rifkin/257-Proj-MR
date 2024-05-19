%% YALMIP OCP relaxation
clc; clear; close all;

p = 3; % x in R^p
n = p + 1; % [x; u] in R^n
x = sdpvar(p, 1); % [s c thdot]
u = sdpvar(1);

kappa = 1; % relaxation order
d = 2*kappa; % degree -- hardcoded monomial methods need to be rewritten if changed
b = monolist([x; u], d); % monomial basis

params.m = 1;
params.l = 1;
params.b = 0.1;
params.g = 9.8;
f = pendulum_dynamics(x, u, params); % pendulum dynamics dx/dt
g = x(1)^2 + x(2)^2 + x(3)^2 + 0.25*u^2; % quadratic cost

S = sdpvar(s(p, d)); % value function gram matrix
J = monolist(x, d)'*S*monolist(x, d); % value function in x

% Constraint for HJB being SOS
HJB = jacobian(J, x) * f + g; % HJB

Q = sdpvar(s(n, d)); % SOS gram matrix
sig = b' * Q * b; % SOS polynomial 

% Constraint for sin^2+cos^2 = 1
h = x(1)^2 + x(2)^2 - 1;
r = sdpvar(s(n, d - degree(h)), 1); 
lam = monolist([x; u], d - degree(h))'*r; % arbitrary polynomial

HJB_constr = coefficients(HJB - sig - lam*h, [x; u]) == 0; % HJB is nonnegative (psatz)

target_x_mon = monomialsN3D2([0, -1, 0]);
target_constr = target_x_mon * S * target_x_mon' == 0; % cost at top of pendulum = 0
 
% Collect constraints
F = [Q>=0, HJB_constr, target_constr];

% set objective (maximize J at x0)
x0_mon = monomialsN3D2([0, 1, 0]); % bottomright position, theta=0
obj = -x0_mon*S*x0_mon';

options = sdpsettings('solver','mosek');
optimize(F,obj,options);

M = dual(F(1)); % moment matrix

regpar = 1e-8;
[V, D] = eig(M);
eigs = diag(D) + regpar; % Tikhonov regularization 

% CD coeff matrix
P = diag(eigs.^(-1/2))*V';
Psq = P'*P; 

%% Offline u computation from CD kernel at sample x
% samples of x
Nxrt = 10;
Nx = Nxrt^p;
x1_vals = linspace(-1, 1, Nxrt);
x2_vals = x1_vals;
x3_vals = x1_vals;

% samples of u
Nu = 100;
u_vals = linspace(-1, 1, Nu);

% create all sample combinations of x and u
[U, X3, X2, X1] = ndgrid(u_vals, x3_vals, x2_vals, x1_vals); % U changes faster than X3 than X2 than X1
X1 = X1(:); X2 = X2(:); X3 = X3(:); U = U(:);

b_eval = monomialsN4D2([X1, X2, X3, U]); % degree 2 monomial vector at all points
q_eval = reshape(sum(b_eval*Psq.*b_eval, 2), Nu, Nx)'; % dot product -> sum to avoid multiplying large matrix
% q_eval is a Nx by Nu matrix where u is compared in each row
[~, idx] = min(q_eval, [], 2); % this will also take outer min as u is sorted in ascending order
u_out = u_vals(idx);


%% Online u compuatation from CD kernel at current X
close all; 
T = 30; %
break_interval = 0.05; % interval to reevaluate u
N_breaks = T/break_interval;

u_traj = [];

x0 = [0, 1, 0];
x_traj = x0;


for i = 1:N_breaks
    % samples of u
    Nu = 1000;
    u_vals = linspace(-1, 1, Nu);
    Nx = 1;

    [U, X3, X2, X1] = ndgrid(u_vals, x0(3), x0(2), x0(1)); 
    X1 = X1(:); X2 = X2(:); X3 = X3(:); U = U(:);
    b_eval = monomialsN4D2([X1, X2, X3, U]); % degree 2 monomial vector at all points
    q_eval = reshape(sum(b_eval*Psq.*b_eval, 2), Nu, Nx)'; % dot product -> sum to avoid multiplying large matrix
    % q_eval is a Nx by Nu matrix where u is compared in each row
    [~, idx] = min(q_eval, [], 2); % this will also take outer min as u is sorted in ascending order
    u_x0 = u_vals(idx);
   
    u_traj = [u_traj, u_x0];

    [~, x_traj_i] = ode45(@(t, x) pendulum_dynamics(x, u_x0, params), [0 break_interval], x0);
    x0 = x_traj_i(end, :);
    x_traj = [x_traj; x0];
end

th_thdot_traj = [atan2(x_traj(:, 1), x_traj(:, 2)), x_traj(:, 3)];
figure;
plot(th_thdot_traj(:, 1), th_thdot_traj(:, 2));
xlabel("$\theta$", 'Interpreter', 'latex');
ylabel("$\dot{\theta}$", 'Interpreter', 'latex');

tt = 0:break_interval:T;
figure;
tiledlayout(3, 1);
ax1 = nexttile;
plot(tt, th_thdot_traj(:, 1));
ylabel("$\theta$", 'Interpreter', 'latex');
grid on;

ax2 = nexttile;
plot(tt, th_thdot_traj(:, 2));
ylabel("$\dot{\theta}$", 'Interpreter', 'latex');
grid on;

ax3 = nexttile;
plot(tt(1:end-1), u_traj);
ylabel("$u$", 'Interpreter', 'latex');
grid on;

