%% YALMIP OCP relaxation
clc; clear; close all;

x = sdpvar(2,1); % [x1, x2]
u = sdpvar(1);

kappa = 0.5; % relaxation order
d = 2*kappa; % degree
b = monolist([x; u], d); % monomial basis

f = [x(2); u]; % double integrator dynamics dx/dt
cost = x(1)^2 + x(2)^2 + u^2; % quadratic cost

S = sdpvar(s(2, d)); % value function gram matrix
J = monolist(x, d)'*S*monolist(x, d); % value function in x

HJB = jacobian(J, x) * f + cost; % HJB

Q0 = sdpvar(s(3, d)); % SOS gram matrix
sig0 = b' * Q0 * b; % SOS polynomial

F = Q0 >=0;
rhs = sig0; 

% State and control inequality constraints 
g = [x(1) + 1;
 -x(1) + 1; 
 x(2) + 1;
-x(2) + 1; 
u + 1;
- u + 1];

for i = 1:length(g)
   deg_sig = floor((d - degree(g(i)))/2);
   Qi = sdpvar(s(3, deg_sig)); 
   F = [F, Qi >= 0];
   sigi = monolist([x; u], deg_sig)'*Qi*monolist([x; u], deg_sig);
   rhs = rhs + sigi*g(i);
end
  
HJB_constr = coefficients(HJB - rhs, [x; u]) == 0; % HJB is nonnegative (psatz)

target_x = [0; 0];
target_constr =  replace(J, x, target_x) == 0; % target constraint 

% Collect constraints
F = [F; HJB_constr, target_constr];

% set objective (maximize J at x0)
x0 = [0.5; 1];
obj = -replace(J, x, x0);

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

% samples of x
Nxrt = 50;
Nx = Nxrt^2;
x1_vals = linspace(-1, 1, Nxrt);
x2_vals = x1_vals;

% samples of u
Nu = 1000;
u_vals = linspace(-1, 1, Nu);

% create all sample combinations of x and u
[U, X2, X1] = ndgrid(u_vals, x2_vals, x1_vals); % U changes faster than X2 than X1
X1 = X1(:); X2 = X2(:); U = U(:);

b_eval = eval_monomials([X1, X2, U], d); % degree d monomial vector at all points
q_eval = reshape(sum(b_eval*Psq.*b_eval, 2), Nu, Nx)'; % dot product -> sum to avoid multiplying large matrix
% q_eval is a Nx by Nu matrix where u is compared in each row
[~, idx] = min(q_eval, [], 2); % this will also take outer min s u is sorted in ascending orderciln
u_out = u_vals(idx);

% calculate correct u from ARE
[X1_noU, X2_noU] = meshgrid(x1_vals, x2_vals); X1_noU = X1_noU(:); X2_noU = X2_noU(:);
u_are_vec = clip(-K_are*[X1_noU, X2_noU]', -1, 1);
u_err = abs(u_are_vec - u_out);

find(u_err > 2/Nu) % this is empty if u_out is closest possible to ARE
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

