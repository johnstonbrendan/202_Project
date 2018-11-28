m = 12.0; %mass of frisbee (grams)
g = 9.81; %gravity m/s^2
CL0 = 1.0; %coefficient of lift (0)
CLa = 1.0; %coefficient of lift (alpha)
CD0 = 1.0; %coefficient of drag (0)
CDa = 1.0; %coefficient of drag (alpha);
alpha = 1.0; %alpha (radians)
alpha_i = 1.0; %alpha initial (radians)
rho = 1.0; %density of fliud (NEED UNITS)
r = 1.0; %radius (m)
y_velocity = 0; % this is a assumption for now (m/s)


syms z(t) x(t) z_velocity(t) x_velocity(t) t X Z Y

ode1 = z_velocity == diff(z);
ode2 = x_velocity == diff(x);
ode3 = x_velocity + lift_force(CL0,CLa,alpha,rho,r,velocity(x_velocity,y_velocity,z_velocity))*z_velocity...
        /velocity(x_velocity,y_velocity,z_velocity) - ...
        drag_force(CD0,CDa,alpha,alpha_i,rho,r,velocity(x_velocity,y_velocity,z_velocity))*x_velocity/...
        velocity(x_velocity,y_velocity,z_velocity) == m*diff(x_velocity);
ode4 = z_velocity + lift_force(CL0,CLa,alpha,rho,r,velocity(x_velocity,y_velocity,z_velocity))*x_velocity...
        /velocity(x_velocity,y_velocity,z_velocity) + ...
        drag_force(CD0,CDa,alpha,alpha_i,rho,r,velocity(x_velocity,y_velocity,z_velocity))*z_velocity/...
        velocity(x_velocity,y_velocity,z_velocity) - m*g == m*diff(z_velocity);
  
[DEVF, Subs] = odeToVectorField(ode1, ode2, ode3, ode4);
ODEfcn = matlabFunction(DEVF, 'Vars', {t, X, Y})
[t, X] = ode45(@(t,Y) ODEfcn(t, X, Y), [0:200],[0;2;0;0])
