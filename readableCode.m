
function solver
    % time step setup
    tstart = 0;
    tend = 20;
    tstep = 0.001;
    n = (tend-tstart)/tstep;
    tspan = linspace(tstart, tend, n);

    % initial conditions
    startHeight = 0; %m
    wind_speed_x = 0;%(m/s)
    wind_speed_z = 0;%(m/s)
    throwV = 15;
    pitch = 10*pi/180;
    xInit = [0;15; 0; 0; startHeight; 0]; %x, vx, y, vy, z, vz, 
    [t, out] = ode45(@(t, out) discODEs(t, out, pitch), tspan, xInit)

    x = out(:,1) +(t*wind_speed_x);
    vx = out(:,2) +(wind_speed_x);
    y = out(:,3);
    vy = out(:,4);
    z = out(:,5) +(t*wind_speed_z);
    vz = out(:,6) +(wind_speed_z);

    showPlots(t, x, y, z, vx, vy, vz, startHeight);
end

function ddt = discODEs(t, out, pitch)
    m = 0.175; %mass of frisbee (kg)
    g = 9.81; %gravity m/s^2
    CL0 = 0.15; %coefficient of lift (0)
    CLa = 1.4; %coefficient of lift (alpha)
    CD0 = 0.08; %coefficient of drag (0)
    CDa = 2.72; %coefficient of drag (alpha);
    alpha_0 = -0.0698; %alpha(diff(x),diff(z)) 0 defined based on physcial aspects of frisbee (radians)
    rho = 1.225; %density of fliud (NEED UNITS)
    r = 0.137; %radius (m)
    roll = 0;
    

    
    % ddt order is [vx, ax, jx, vz, az, jz] where jx, jz are jerk of x
    x = out(1);
    vx = out(2);
    y = out(3);
    vy = out(4);
    z = out(5);
    vz = out(6);
    
    % frisbee position and velocity vector wrt XYZ
    p = [x y z];
    v = [vx vy vz]; 
    
    
    lift_vect = cross(v, [0 1 tan(roll)]);
    lift_uvect = lift_vect/norm(lift_vect);
    
    drag_vect = -1*v;
    drag_uvect = drag_vect/norm(drag_vect);
    
    lift_force = calc_lift_force(CL0, CLa, alpha(vx, vz, pitch), rho, r, velocity(vx, vy, vz));
    lift = lift_force*lift_uvect;
    
    drag_force = calc_drag_force(CD0, CDa, alpha(vx, vz, pitch), alpha_0, rho, r, velocity(vx, vy, vz));
    drag = drag_force*drag_uvect;
    
    liftX = lift(1);
    liftY = lift(2);
    liftZ = lift(3);
    
    dragX = drag(1);
    dragY = drag(2);
    dragZ = drag(3);
    % X = [1 0 0], Y = [0 1 0], Z = [0 0 1]
    
    ddt = zeros(size(out));
    ddt(1) = vx;
    ddt(2) = (liftX + dragX)/m;
    ddt(3) = vy;
    ddt(4) = (liftY + dragY)/m;
    ddt(5) = vz;
    ddt(6) = (liftZ + dragZ - m*g)/m;
end

function showPlots(t, x, y, z, vx, vy, vz, startHeight)
    figure('Name','x y z plot')
    plot3(x, y,z)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    xlim([0 20])
    ylim([0 20])
    zlim([0 20])
    figure('Name','z t plot')
    plot(t,z)
    xlabel('t')
    ylabel('z')
    ylim ([0 5])
    
    figure('Name','x t plot')
    plot(t,x)
    xlabel('t')
    ylabel('x')
    
     
    figure('Name','y t plot')
    plot(t,y)
    xlabel('t')
    ylabel('y')
end
    
