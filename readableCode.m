function solver
    %%%%%%%%% initial conditions %%%%%%%%
    startHeight = 2; %m
    wind_speed = [0 0 0]; %x,y,z windspeed (m/s)
    throwV = 14; % magnitude of the throw velocity (m/s)
    initPitch =10*pi/180; % pitch of initial throw (rad)
    initRoll = 0*pi/180; %roll of initial throw (rad)
    spinRate =37; % spin throw rate of frisbee (rad/s)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % time step setup
    tstart = 0; %seconds
    tend = 200; %this can be large as the calculation will stop when z = 0
    tstep = 0.0001;
    n = (tend-tstart)/tstep;
    tspan = linspace(tstart, tend, n);

    
    xInit = [0; throwV*cos(initPitch); 0; 0; startHeight; throwV*sin(initPitch); initRoll; initPitch]; %x, vx, y, vy, z, vz, r, p
    Opt = odeset('Events', @detectGround);
    [t, out] = ode45(@(t, out) discODEs(t, out, wind_speed, spinRate), tspan, xInit, Opt);

    x = out(:,1) +(t*wind_speed(1));
    vx = out(:,2) + wind_speed(1);
    y = out(:,3) + (t*wind_speed(2));
    vy = out(:,4) + wind_speed(2);
    z = out(:,5) +(t*wind_speed(3));
    vz = out(:,6) + wind_speed(3);
    r = out(:, 7).*180/pi; % roll state (degrees)
    p = out(:, 8).*180/pi; % pitch state (degrees)

   showPlots(t, x, y, z, vx, vy, vz, r, p, startHeight);
   
end

function ddt = discODEs(t, out, wind_speed, spinRate)
    m = 0.175; %mass of frisbee (kg)
    g = 9.81; %gravity (m/s^2)
    CL0 = 0.15; %coefficient of lift (0)
    CLa = 1.4; %coefficient of lift (alpha)
    CD0 = 0.08; %coefficient of drag (0)
    CDa = 2.72; %coefficient of drag (alpha);
    alpha_0 = -0.0698; % 0 defined based on physcial aspects of frisbee (radians)
    rho = 1.225; %density of fliud (kg/m^3)
    rd = 0.137; %radius (m)
    Iz = 1/8*m*(rd*2)^2; % kg(m^2)
    

    
    % ddt order is [vx, ax, vy, ay, vz, az, rdot, pdot] 
    x = out(1);
    vx = out(2);
    y = out(3);
    vy = out(4);
    z = out(5);
    vz = out(6);
    r = out(7);
    p = out(8);
    
    % frisbee position and velocity vector wrt XYZ
    pos = [x y z];
    v = [vx vy vz]; 
    v = v - wind_speed; % find free stream velocity, v frisbee rel to air
    
    % global to body transformation matrix
    globalToBody = transpose([1 0 0;0 cos(r) -sin(r);0 sin(r) cos(r)]*[cos(p) 0 sin(p); 0 1 0; -sin(p) 0 cos(p)]);
    
    % transform global axes to body axes
    bodyx = (globalToBody*transpose([1 0 0])).';
    bodyy = (globalToBody*transpose([0 1 0])).';
    bodyz = (globalToBody*transpose([0 0 1])).';
    
    %Calculate lift and drag unit vector
    lift_vect = cross(v, bodyy);
    lift_uvect = lift_vect/norm(lift_vect);
    drag_vect = -1*v;
    drag_uvect = drag_vect/norm(drag_vect);
    
    %Calculate lift and drag force
    lift_force = calc_lift_force(CL0, CLa, alpha(vx, vz, p), rho, rd, velocity(vx, vy, vz));
    lift = lift_force*lift_uvect;
    drag_force = calc_drag_force(CD0, CDa, alpha(vx, vz, p), alpha_0, rho, rd, velocity(vx, vy, vz));
    drag = drag_force*drag_uvect;
    
    %Get lift and drag force components
    liftX = lift(1);
    liftY = lift(2);
    liftZ = lift(3);
    liftBodyz = (dot(lift,bodyz)/(norm(bodyz))^2)*bodyz;
    
    dragX = drag(1);
    dragY = drag(2);
    dragZ = drag(3);
    dragBodyz = (dot(drag,bodyz)/norm(bodyz)^2)*bodyz;
    % X = [1 0 0], Y = [0 1 0], Z = [0 0 1]
    
    % PITCHING MOMENT SECTION % 
    
    % assume clockwise rotation
    totalBodyzForce = norm(liftBodyz + dragBodyz);
    
    % project v vector onto plane
    acDirVector = (v - (dot(v, bodyz)/(norm(bodyz))^2)*bodyz);
    acPosVector = 0.12*rd*2*acDirVector/norm(acDirVector);
    
    % use scalar dot product formula to find angle to bodyx axis
    nu = acos(dot(bodyx, acPosVector))/(norm(bodyx)*norm(acPosVector));
    
    % find components of the position vector to ac on bodyx, bodyy axis
    acPosVectorx = acPosVector*cos(nu);
    acPosVectory = acPosVector*sin(nu);
    
    % find negative pitching moment by using sin(angle) as distance
    pMom = -norm((totalBodyzForce * acPosVectory)/globalToBody);
    
    % find rolling moment by using cos(angle) as distance
    rMom = norm((totalBodyzForce * acPosVectorx)/globalToBody);
    
    ddt = zeros(size(out));
    ddt(1) = vx;
    ddt(2) = (liftX + dragX)/m;
    ddt(3) = vy;
    ddt(4) = (liftY + dragY)/m;
    ddt(5) = vz;
    ddt(6) = (liftZ + dragZ - m*g)/m;
    ddt(7) = 2*rMom/(Iz * spinRate);
    ddt(8) = 2*pMom/(Iz * spinRate);
end

function showPlots(t, x, y, z, vx, vy, vz, r, p, startHeight)
    close all;
    %x,y,z plot
    figure('Name','x y z plot');
    plot3(x, y,z);
    x_limit = xlim();
    y_limit = ylim();
    z_limit = zlim();
    max_lim = max([z_limit(2) y_limit(2) x_limit(2)]);
    axis([0 max_lim 0 max_lim 0 max_lim]);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');
    %z vs time plot 
    figure('Name','z t plot');
    plot(t,z);
    xlabel('t (s)');
    ylabel('z (m)');
    %x vs time plot
    figure('Name','x t plot');
    plot(t,x);
    xlabel('t (s)');
    ylabel('x (m)');
    %y vs time plot
    figure('Name','y t plot');
    plot(t,y);
    xlabel('t (s)');
    ylabel('y (m)');
    %roll vs time plot
    figure('Name', 'r t plot');
    plot(t, r);
    xlabel('t (s)');
    ylabel('roll (degrees)');
    %pitch vs time plot
    figure('Name', 'p t plot');
    plot(t, p);
    xlabel('t (s)');
    ylabel('pitch (degrees)');
end
    
function [value, isterminal, direction] = detectGround(t, out)
    value = (out(5) < 0);
    isterminal = 1;
    direction = 0;
end
