function f = pendulum_dynamics(x, u, params)
    f = [x(2)*x(3); 
        -x(1)*x(3); 
        -(params.b*x(3) + params.m*params.g*params.l*x(1))/(params.m*params.l^2)];
    
    f = f + [0; 0; 1/(params.m*params.l^2)]*u;
end