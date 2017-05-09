function [ tout, xout ] = integrator_rk4(dt,dx,t,x)
    
    % function xnext = integrator_rk4(dt,dx,t,x)
    %
    % This function performs a single Runge-Kutta update step
    %
    % input: 
    % dt: timestep
    % dx: function handle to dx(t,x)
    % t: current time
    % x: current state (column vector)
    %
    % output:
    % tout: t + dt
    % xout: state at time t + dt
    %
    tout = t + dt;
    dt_half = 0.5*dt;
    
    k1 = dx(t,x);
    k2 = dx(t+dt_half,x+dt_half*k1);
    k3 = dx(t+dt_half,x+dt_half*k2);
    k4 = dx(tout,x+dt*k3);
    xout = x + dt*(k1+2*k2+2*k3+k4)/6;
end