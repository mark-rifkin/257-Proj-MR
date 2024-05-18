%% YALMIP OCP relaxation
clc; clear; close all;

x = sdpvar(2,1); % [x1, x2]
u = sdpvar(1);

d = 1; % relaxation order
b = monolist([x; u], d); % monomial basis


f = [x(2); u]; % double integrator dynamics dx/dt
g = x(1)^2 + x(2)^2 + u^2; % quadratic cost

S = sdpvar(s(2, d)); % value function gram matrix
J = monolist(x, d)'*S*monolist(x, d); % value function in x

HJB = jacobian(J, x) * f + g; % HJB

Q = sdpvar(s(3, d)); % SOS gram matrix
sig = b' * Q * b; % SOS polynomial equal to HJB

target = [1; 0; 0;]' * S * [1; 0; 0]; % target constraint (wasy to use subs?)
 
% Collect constraints
F = [Q>=0, coefficients(HJB - sig, [x; u]) == 0, target == 0];

% set objective (maximize J at x0)
x0 = [1; 0.5; 1];
obj = -x0'*S*x0;

options = sdpsettings('solver','mosek');
optimize(F,obj,options);

S_val = value(S);

% check S from ARE
[K_are, S_are] = lqr([0, 1; 0, 0], [0; 1], [1, 0; 0, 1], 1);

M = dual(F(1)); % moment matrix

%% U computation from CD kernel
regpar = 1e-8;
[V, D] = eig(M);
eigs = diag(D) + regpar; % Tikhonov regularization 

% CD coeff matrix
P = diag(eigs.^(-1/2))*V';
Psq = P'*P; 

x = sym('x', [2 1]);
u = sym('u');
b = monomials([x; u], 0:d);

% samples of x
Nxrt = 100;
Nx = Nxrt^2;
x1_vals = linspace(-1, 1, Nxrt);
x2_vals = x1_vals;

% samples of u
Nu = 1000;
u_vals = linspace(-1, 1, Nu);

% create all sample combinations of x and u
[X2, U, X1] = meshgrid(x2_vals, u_vals, x1_vals); % U changes faster than X2 than X1
X1 = X1(:); X2 = X2(:); U = U(:);

b_eval = [ones(Nx*Nu, 1), X1, X2, U]; % 1-degree monomial vector at all points
q_eval = reshape(sum(b_eval*Psq.*b_eval, 2), Nu, Nx)'; % dot product -> sum to avoid multiplying large matrix
% q_eval is a Nx by Nu matrix where u is compared in each row
[~, idx] = min(q_eval, [], 2); % this will also take outer min s u is sorted in ascending orderc
u_out = u_vals(idx);

% calculate correct u from ARE
[X1_noU, X2_noU] = meshgrid(x1_vals, x2_vals); X1_noU = X1_noU(:); X2_noU = X2_noU(:);
u_are_vec = clip(-K_are*[X1_noU, X2_noU]', -1, 1);
u_err = abs(u_are_vec - u_out);

find(u_err > 2/Nu) % this is empty so u_out is closest possible to ARE
%% SOSTools (to check YALMIP)
clc; clear; close all;
var = mpvar('var', [3 1]);
x = var(1:2);
u = var(end);

prog = sosprogram(var);

degree = 2;
[prog, J] = sospolyvar(prog, monomials(x, 0:degree));

f = [x(2); u]; % double integrator dynamics dx/dt
g = x(1)^2 + x(2)^2 + u^2; % quadratic cost

% HJB inequality constraint
HJB = diff(J, x)*f + g;
prog = sosineq(prog, HJB);

% Target value 0 equality constraint
J_0 = subs(J, x, [0; 0]);
prog = soseq(prog, J_0);

% set objective (maximize J at x0)
x0 = [0.5; 1];
J_x0 = subs(J, x, x0);
prog = sossetobj(prog, -J_x0);

% solve for J
options.solver = 'mosek';
prog = sossolve(prog, options);
Jstar = sosgetsol(prog, J);

