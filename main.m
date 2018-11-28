m = 12.0; %mass of frisbee (grams)
g = 9.81; %gravity m/s^2
CL0 = 1.0; %coefficient of lift (0)
CLa = 1.0; %coefficient of lift (alpha)
CD0 = 1.0; %coefficient of drag (0)
CDa = 1.0; %coefficient of drag (alpha);
a = 1.0; %alpha (radians)
rho = 1.0; %density of fliud (NEED UNITS)
r = 1.0; %radius (m)
y_velocity = 0; % this is a assumption for now (m/s)


syms z(t),x(t),z_velocity(t),x_velocity(t)
ode1 = diff(z_velocity) == z;
ode2 = diff(x_veloctiy) == x;
ode3 = x_velocity + lift_force(CL0,CLa,alpha,rho,r,v)*z_velocity...
        /velocity(x_velocity,y_velocity,z_velocity) - ...
        lift_drag(CD0,CDa,alpha,alpha_i,rho,r,v)*x_velocity/...
        velocity(x_velocity,y_velocity,z_velocity) == m*diff(x_velocity);
ode4 = z_velocity + lift_force(CL0,CLa,alpha,rho,r,v)*x_velocity...
        /velocity(x_velocity,y_velocity,z_velocity) + ...
        lift_drag(CD0,CDa,alpha,alpha_i,rho,r,v)*z_velocity/...
        velocity(x_velocity,y_velocity,z_velocity) - m*g == m*diff(x_velocity);