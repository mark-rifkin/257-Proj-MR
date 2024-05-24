function [f, f2] = pendulum_dynamics(x, u, params)
    thdot = 4*pi*x(3);
    f2 = [0; 0; 1/(params.m*params.l^2)];
    f = [x(2)*thdot; 
        -x(1)*thdot; 
        -(params.b*thdot + params.m*params.g*params.l*x(1))/(params.m*params.l^2)]...
    + f2*u;
    
end