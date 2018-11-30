
function solver
    % time step setup
    tstart = 0;
    tend = 100;
    tstep = 0.001;
    n = (tend-tstart)/tstep;
    tspan = linspace(tstart, tend, n);

    % initial conditions
    startHeight = 2; %m
    wind_speed = [0 0 0]; %x,y,z windspeed (m/s)
    throwV = 15;
    initPitch = 10*pi/180;
    spinRate = 60*pi;
    xInit = [0; throwV*cos(initPitch); 0; 0; startHeight; throwV*sin(initPitch); 0*pi/180; initPitch]; %x, vx, y, vy, z, vz, r, p
    Opt = odeset('Events', @detectGround);
    [t, out] = ode45(@(t, out) discODEs(t, out, wind_speed, spinRate), tspan, xInit, Opt)

    x = out(:,1) +(t*wind_speed(1));
    vx = out(:,2) + wind_speed(1);
    y = out(:,3) + (t*wind_speed(2));
    vy = out(:,4) + wind_speed(2);
    z = out(:,5) +(t*wind_speed(3));
    vz = out(:,6) + wind_speed(3);
    r = out(:, 7).*180/pi; % roll state (degrees)
    p = out(:, 8).*180/pi; % pitch state (degrees)

   showPlots(t, x, y, z, vx, vy, vz, r, p, startHeight);
   
   %ylim ([0 startHeight+4]);
   %xlim ([0 10]);
   
end

function ddt = discODEs(t, out, wind_speed, spinRate)
    m = 0.175; %mass of frisbee (kg)
    g = 9.81; %gravity m/s^2
    CL0 = 0.15; %coefficient of lift (0)
    CLa = 1.4; %coefficient of lift (alpha)
    CD0 = 0.08; %coefficient of drag (0)
    CDa = 2.72; %coefficient of drag (alpha);
    alpha_0 = -0.0698; %alpha(diff(x),diff(z)) 0 defined based on physcial aspects of frisbee (radians)
    rho = 1.225; %density of fliud (NEED UNITS)
    rd = 0.137; %radius (m)
    Iz = 1/8*m*(rd*2)^2;
    

    
    % ddt order is [vx, ax, jx, vz, az, jz] where jx, jz are jerk of x
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
    
    globalToBody = transpose([1 0 0;0 cos(r) -sin(r);0 sin(r) cos(r)]*[cos(p) 0 sin(p); 0 1 0; -sin(p) 0 cos(p)]);
    bodyx = (globalToBody*transpose([1 0 0])).';
    bodyy = (globalToBody*transpose([0 1 0])).';
    bodyz = (globalToBody*transpose([0 0 1])).';
    
    lift_vect = cross(v, bodyy);
    lift_uvect = lift_vect/norm(lift_vect);
    
    drag_vect = -1*v;
    drag_uvect = drag_vect/norm(drag_vect);
    
    lift_force = calc_lift_force(CL0, CLa, alpha(vx, vz, p), rho, rd, velocity(vx, vy, vz));
    lift = lift_force*lift_uvect;
    
    drag_force = calc_drag_force(CD0, CDa, alpha(vx, vz, p), alpha_0, rho, rd, velocity(vx, vy, vz));
    drag = drag_force*drag_uvect;
    
    liftX = lift(1);
    liftY = lift(2);
    liftZ = lift(3);
    liftBodyz = (dot(lift,bodyz)/(norm(bodyz))^2)*bodyz;
    
    dragX = drag(1);
    dragY = drag(2);
    dragZ = drag(3);
    dragBodyz = (dot(drag,bodyz)/norm(bodyz)^2)*bodyz;
    % X = [1 0 0], Y = [0 1 0], Z = [0 0 1]
    
    % assume clockwise rotation
    totalBodyzForce = norm(liftBodyz + dragBodyz);
    
    % project v vector onto plane
    acDirVector = (v - (dot(v, bodyz)/(norm(bodyz))^2)*bodyz)
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
    close all
    figure('Name','x y z plot')
    plot3(x, y,z)
    xlabel('x')
    ylabel('y')
    zlabel('z')

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
    
    figure('Name', 'r t plot')
    plot(t, r)
    xlabel('t')
    ylabel('r')
    
    figure('Name', 'p t plot')
    plot(t, p)
    xlabel('t')
    ylabel('p')
    
    %title = ['Windspeed is ',num2str(wind_speed(3)),' in z'];
   %figure('Name',title)
   %plot (t,z);
   %xlabel ('t');
   %ylabel ('z');
end
    
function [value, isterminal, direction] = detectGround(t, out)
    value = (out(5) < 0);
    isterminal = 1;
    direction = 0;
end
    
