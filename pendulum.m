%% YALMIP OCP relaxation
clc; clear; close all;
addpath(genpath(pwd));

p = 3; % x in R^p
n = p + 1; % [x; u] in R^n
x = sdpvar(p, 1); % [s c thdot]
u = sdpvar(1);

kappa = 1; % relaxation order
d = 2*kappa; % degree
b = monolist([x; u], d); % monomial basis

params.m = 1;
params.l = 1;
params.b = 0.1;
params.g = 9.8;
[f, f2] = pendulum_dynamics(x, u, params); % pendulum dynamics dx/dt
rho = 0.01;
% target_x = [0; -1; 0];
target_x = [0; 0; 0];

% x_err = x - target_x;
cost = x'*eye(p)*x + u^2; % quadratic cost

S = sdpvar(s(p, d)); % value function gram matrix
J = monolist(x, d)'*S*monolist(x, d); % value function in x

HJB = -rho*J + jacobian(J, x) * f + cost; % HJB

% State and control inequality constraints 

g = [1- [x; u]'*[x; u]];

Q = cell(1+length(g), 1);
sig = cell(1+length(g), 1);

% g0 constraint
Q{1} = sdpvar(s(n, d)); % SOS gram matrix
sig{1} = b' * Q{1} * b; % SOS polynomial 

F = Q{1}>= 0;
rhs = sig{1};

for i = 1:length(g)
    % generate SOS multiplier
   deg_sig = floor((d - degree(g(i)))/2);
   Q{i+1} = sdpvar(s(n, deg_sig)); 
   sig{i+1} = monolist([x; u], deg_sig)'*Q{i+1}*monolist([x; u], deg_sig);

   F = [F, Q{i+1} >=0];  % sig(x) is SOS
  
   rhs = rhs + sig{i+1}*g(i);
end

% State equality constraint
h = x(1)^2 + x(2)^2 - 1;
deg_lam = d - degree(h);
r = sdpvar(s(n, deg_lam), 1);
lam = monolist([x; u], deg_lam)'*r;
rhs = rhs + lam*h;

HJB_constr = coefficients(HJB - rhs, [x; u]) == 0; % HJB is nonnegative (psatz)

target_constr = replace(J, x, target_x) == 0; % cost at top of pendulum = 0
 
% Collect constraints
F = [F, HJB_constr, target_constr];

% set objective (maximize J at x0)
x0 = [sin(0.1*pi); 1+cos(0.1*pi); 0];
obj = -replace(J, x, x0);

options = sdpsettings('solver','mosek');
optimize(F,obj,options);

S_opt = value(S);

M = dual(F(1)); % moment matrix

regpar = 1e-8;
[V, D] = eig(M);
eigs = diag(D) + regpar; % Tikhonov regularization 

% CD coeff matrix
P = diag(eigs.^(-1/2))*V';
Psq = P'*P; 

%% Estimate u from CD kernel of moment matrix at sample x

N_th = 5;
N_thdot = 3;
Nx = N_th*N_thdot;

th_vals = linspace(-pi, pi, N_th);
thdot_vals = linspace(-0.5*pi, 0.5*pi, N_thdot);

% samples of u
Nu = 101;
u_vals = linspace(-1, 1, Nu);

[U, Th, Thdot] = ndgrid(u_vals, th_vals, thdot_vals); % all sample combinations of x and u
U = U(:);
X1 = sin(Th(:)); 
X2 = cos(Th(:))-1;
X3 = Thdot(:);

b_eval = eval_monomials([X1, X2, X3, U], d); % degree d monomial vector at all points
q_eval = reshape(sum(b_eval*Psq.*b_eval, 2), Nu, Nx)'; % dot product -> sum to avoid multiplying large matrix
% q_eval is a Nx by Nu matrix where u is compared in each row
[~, idx] = min(q_eval, [], 2); % this will also take outer min as u is sorted in ascending order
u_out = u_vals(idx);

%% Compute u directly from dJ/dx at sample x
[Th, Thdot] = ndgrid(th_vals, thdot_vals);
X1 = sin(Th(:)); 
X2 = cos(Th(:))-1;
X3 = Thdot(:);

mu = -0.5*f2'*jacobian(J, x)';

u_direct = zeros(size(X3));

for i = 1:length(u_direct)
    u_direct(i) = clip(replace(mu, x, [X1(i); X2(i); X3(i)]), -1, 1);
end

find(u_out-u_direct' > 2/Nu)

%% Optimize at different x0 to see if different J generated
th_vals = [-pi, -3.141, -3.14, -3.1];
thdot_vals = 0;
[Th, Thdot] = ndgrid(th_vals, thdot_vals);
X1 = sin(Th(:)); 
X2 = cos(Th(:))-1;
X3 = Thdot(:);

x0 = [X1(:), X2(:), X3(:)];
S_mat = zeros(s(p, d)^2, size(x0, 1)); % stores all S_opt matrices
for i = 1:size(x0, 1)
    obj = -replace(J, x, x0(i, :));
    optimize(F,obj,options);
    S_mat(:, i) = reshape(value(S), s(p, d)^2, 1);
end

% compare these programmatically? 
%% Online u compuatation from CD kernel at current X
close all; 
T = 30; %
break_interval = 0.05; % interval to reevaluate u
N_breaks = T/break_interval;

u_traj = [];

x0 = [0, 1, 0.1];
x_traj = x0;


% samples of u
Nu = 101;
u_vals = linspace(-1, 1, Nu);
Nx = 1;

for i = 1:N_breaks

    [U, X3, X2, X1] = ndgrid(u_vals, x0(3), x0(2), x0(1)); 
    X1 = X1(:); X2 = X2(:); X3 = X3(:); U = U(:);
    b_eval = eval_monomials([X1, X2, X3, U], d); % degree d monomial vector at all points
    q_eval = reshape(sum(b_eval*Psq.*b_eval, 2), Nu, Nx)'; % dot product -> sum to avoid multiplying large matrix
    % q_eval is a Nx by Nu matrix where u is compared in each row
    [~, idx] = min(q_eval, [], 2); % this will also take outer min as u is sorted in ascending order
    u_x0 = u_vals(idx);
   
    u_traj = [u_traj, u_x0];

    [~, x_traj_i] = ode45(@(t, x) pendulum_dynamics(x, u_x0, params), [0 break_interval], x0);
    x0 = x_traj_i(end, :);
    x_traj = [x_traj; x0];
end

th_thdot_traj = [atan2(x_traj(:, 1), x_traj(:, 2)-1), x_traj(:, 3)];
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

