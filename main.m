m = 12.0; %mass of frisbee (grams)
g = 9.81; %gravity m/s^2
CL0 = 1.0; %coefficient of lift (0)
CLa = 1.0; %coefficient of lift (alpha)
CD0 = 1.0; %coefficient of drag (0)
CDa = 1.0; %coefficient of drag (alpha);
alpha_i = 1.0; %alpha(x_velocity,z_velocity) initial (radians)
rho = 1.0; %density of fliud (NEED UNITS)
r = 1.0; %radius (m)
y_velocity = 0; % this is a assumption for now (m/s)


syms z(t) x(t) z_velocity(t) x_velocity(t)
ode1 = diff(z) == z_velocity;
ode2 = diff(x) == x_velocity;
ode3 = diff(x_velocity) == (x_velocity + lift_force(CL0,CLa,alpha(x_velocity,z_velocity),rho,r,velocity(x_velocity,y_velocity,z_velocity))*z_velocity...
        /velocity(x_velocity,y_velocity,z_velocity) - ...
        drag_force(CD0,CDa,alpha(x_velocity,z_velocity),alpha_i,rho,r,velocity(x_velocity,y_velocity,z_velocity))*x_velocity/...
        velocity(x_velocity,y_velocity,z_velocity))/m;
ode4 = diff(z_velocity) == (z_velocity + lift_force(CL0,CLa,alpha(x_velocity,z_velocity),rho,r,velocity(x_velocity,y_velocity,z_velocity))*x_velocity...
        /velocity(x_velocity,y_velocity,z_velocity) + ...
        drag_force(CD0,CDa,alpha(x_velocity,z_velocity),alpha_i,rho,r,velocity(x_velocity,y_velocity,z_velocity))*z_velocity/...
        velocity(x_velocity,y_velocity,z_velocity) - m*g)/m;
    
cond1 = x_velocity(0) == 0;
cond2 = z_velocity(0) == 2;

    
odes = [ode1;ode3;ode4];
conds = [cond1; cond2];
S = dsolve(odes,conds);
[xSol(t),zSol(t)] = dsolve(odes,conds);
%zSol(t) = S.x_velocity;
%zSol(t) = S.z_velocity;
clf
fplot(xSol)
hold on
fplot(zSol)
grid on
legend('xSol','zSol','Location','best')
