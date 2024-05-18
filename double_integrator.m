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

% generate test points x
npoints = 10;
[X1, X2] = meshgrid(linspace(-1, 1, npoints), linspace(-1, 1, npoints));
X = [X1(:), X2(:)];
Nx = size(X, 1);

% samples of u
Nu = 101;
u_vals = linspace(-1, 1, Nu);

% select u_out from samples of u for test x 
u_out = zeros(1, Nx);
for i = 1:Nx
    b_eval = double(subs(b, {x(1), x(2), u}, {X(i, 1), X(i, 2), u_vals})); % evaluate monomial basis at selected x and all u

    q_eval = diag(b_eval.'*Psq*b_eval); % 1 by Nu -- evaluates q(x, u) for each u and given x
    min_idx = find(q_eval == min(q_eval)); % indices of min(q(x, u))
    u_out(i) = min(u_vals(min_idx)); % min(argmin u q(x,u)) -- outer minimization breaks ties
    
end

u_are_vec = clip(-K_are*X', -1, 1);
u_err = abs(u_are_vec - u_out);
u_err > 2/Nu % this is all 0 so u_out is always closest possible to u from ARE

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

