clc; clear; close all;
%% YALMIP
x = sdpvar(2,1); % [x1, x2]
u = sdpvar(1);

f = [x(2); u]; % double integrator dynamics dx/dt
g = x(1)^2 + x(2)^2 + u^2; % quadratic cost

S = sdpvar(3,3); % value function in terms of x
J = monolist(x, 1)'*S*monolist(x, 1);

HJB = jacobian(J, x) * f + g; % HJB

Q = sdpvar(4, 4);
sig = monolist([x; u], 1)' * Q * monolist([x; u], 1);

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
[~, S_are] = lqr([0, 1; 0, 0], [0; 1], [1, 0; 0, 1], 1);


%% SOSTools (to check)
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

