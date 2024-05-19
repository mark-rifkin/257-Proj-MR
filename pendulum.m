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
g = x(1)^2 + x(2)^2 + x(3)^2 + u^2; % quadratic cost

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

target_x_mon = monomialsN3D2([0, 1, 0]);
target_constr = target_x_mon * S * target_x_mon' == 0; % cost at top of pendulum = 0
 
% Collect constraints
F = [Q>=0, HJB_constr, target_constr];

% set objective (maximize J at x0)
x0_mon = monomialsN3D2([0, -1, 0]);
obj = -x0_mon*S*x0_mon';

options = sdpsettings('solver','mosek');
optimize(F,obj,options);

M = dual(F(1)); % moment matrix

%% U computation from CD kernel
regpar = 1e-8;
[V, D] = eig(M);
eigs = diag(D) + regpar; % Tikhonov regularization 

% CD coeff matrix
P = diag(eigs.^(-1/2))*V';
Psq = P'*P; 

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

%% SOSTools (to check YALMIP gets same objective)
clc; clear; close all;
var = mpvar('var', [4 1]);
x = var(1:3);
u = var(end);

prog = sosprogram(var);

degree = 2;
[prog, J] = sospolyvar(prog, monomials(x, 0:degree));

params.m = 1;
params.l = 1;
params.b = 0.1;
params.g = 9.8;

f = pendulum_dynamics(x, u, params); % double integrator dynamics dx/dt
g = x(1)^2 + x(2)^2 + x(3)^2 + u^2; % quadratic cost

% HJB inequality constraint
HJB = diff(J, x)*f + g;
prog = sosineq(prog, HJB);

% trig equality constraint
h = x(1)^2 + x(2)^2 - 1;
deg_h = h.maxdeg;
kappa = 2;
[prog,lam] = sospolyvar(prog,monomials(var,0:(2*kappa-deg_h)));
prog = soseq(prog, lam*h);

% Target value 0 equality constraint
J_0 = subs(J, x, [0; 1; 0]);
prog = soseq(prog, J_0);

% set objective (maximize J at x0)
x0 = [0; -1; 0];
J_x0 = subs(J, x, x0);
prog = sossetobj(prog, -J_x0);

% solve for J
options.solver = 'mosek';
prog = sossolve(prog, options);
Jstar = sosgetsol(prog, J);

