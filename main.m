m = 12.0; %mass of frisbee (grams)
g = 9.81; %gravity m/s^2
CL0 = 1.0; %coefficient of lift (0)
CLa = 1.0; %coefficient of lift (alpha)
CD0 = 1.0; %coefficient of drag (0)
CDa = 1.0; %coefficient of drag (alpha);
alpha_i = 1.0; %alpha(diff(x),diff(z)) initial (radians)
rho = 1.0; %density of fliud (NEED UNITS)
r = 1.0; %radius (m)
y_velocity = 0; % this is a assumption for now (m/s)


syms z(t) x(t) t Y
ode1 = diff(x, 2) == (diff(x) + lift_force(CL0,CLa,alpha(diff(x),diff(z)),rho,r,velocity(diff(x),y_velocity,diff(z)))*diff(z)...
        /velocity(diff(x),y_velocity,diff(z)) - ...
        drag_force(CD0,CDa,alpha(diff(x),diff(z)),alpha_i,rho,r,velocity(diff(x),y_velocity,diff(z)))*diff(x)/...
        velocity(diff(x),y_velocity,diff(z)))/m;
ode2 = diff(z, 2) == (diff(z) + lift_force(CL0,CLa,alpha(diff(x),diff(z)),rho,r,velocity(diff(x),y_velocity,diff(z)))*diff(x)...
        /velocity(diff(x),y_velocity,diff(z)) + ...
        drag_force(CD0,CDa,alpha(diff(x),diff(z)),alpha_i,rho,r,velocity(diff(x),y_velocity,diff(z)))*diff(z)/...
        velocity(diff(x),y_velocity,diff(z)) - m*g)/m;
  
[DEVF, Subs] = odeToVectorField(ode1, ode2);
ODEfcn = matlabFunction(DEVF, 'Vars', {t, Y});
[t, X] = ode45(@(t,Y) ODEfcn(t, Y), [0 200],[0;0;0;4]);

plot(t, X)


