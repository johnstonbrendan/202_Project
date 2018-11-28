m = 0.175; %mass of frisbee (kg)
g = 9.81; %gravity m/s^2
CL0 = 0.15; %coefficient of lift (0)
CLa = 1.4; %coefficient of lift (alpha)
CD0 = 0.08; %coefficient of drag (0)
CDa = 2.72; %coefficient of drag (alpha);
alpha_i = 0; %alpha(diff(x),diff(z)) initial (radians)
rho = 1.225; %density of fliud (NEED UNITS)
r = 0.137; %radius (m)
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
[t, X] = ode45(@(t,Y) ODEfcn(t, Y), [0 100],[0;0;0;4]); %z, dZ, x, dX

plot(t, X(:,1))



