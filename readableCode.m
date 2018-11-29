
function solver
    startHeight = 2;
    % time step setup
    tstart = 0;
    tend = 4;
    tstep = 0.01;
    n = (tend-tstart)/tstep;
    tspan = linspace(tstart, tend, n);

    % initial conditions
    xInit = [0; 3; startHeight; 0]; %x, vx, z, vz, 
    [t, out] = ode45(@discODEs, tspan, xInit)

    x = out(:,1);
    vx = out(:,2);
    z = out(:,3);
    vz = out(:,4);

    figure('Name','x z t plot')
    plot3(t,x,z)
    xlabel('t')
    ylabel('x')
    zlabel('z')
    zlim([0 startHeight])
    figure('Name','z t plot')
    plot(t,z)
    xlabel('t')
    ylabel('z')
    ylim ([0 startHeight])
    
    figure('Name','x t plot')
    plot(t,x)
    xlabel('t')
    ylabel('x')
    
    %plot stuff here
end

function ddt = discODEs(t, out)
    m = 0.175; %mass of frisbee (kg)
    g = 9.81; %gravity m/s^2
    CL0 = 0.15; %coefficient of lift (0)
    CLa = 1.4; %coefficient of lift (alpha)
    CD0 = 0.08; %coefficient of drag (0)
    CDa = 2.72; %coefficient of drag (alpha);
    alpha_i = 0; %alpha(diff(x),diff(z)) initial (radians)
    rho = 1.225; %density of fliud (NEED UNITS)
    r = 0.137; %radius (m)
    pitch = 0;
    y_velocity = 0; % this is a assumption for now (m/s)
    
    % ddt order is [vx, ax, jx, vz, az, jz] where jx, jz are jerk of x
    x = out(1);
    vx = out(2);
    z = out(3);
    vz = out(4);
    
    
    ddt = zeros(size(out));
    ddt(1) = vx;
    ddt(2) = ((lift_force(CL0, CLa, alpha(vx, vz, pitch), rho, r, velocity(vx, y_velocity, vz))*-vz)/velocity(vx, y_velocity, vz)...
        - drag_force(CD0, CDa, alpha(vx, vz, pitch), alpha_i, rho, r, velocity(vx, y_velocity, vz))*vx/velocity(vx, y_velocity, vz))/m;
    ddt(3) = vz;
    ddt(4) = (lift_force(CL0, CLa, alpha(vx, vz, pitch), rho, r, velocity(vx, y_velocity, vz))*vx/velocity(vx, y_velocity, vz) ...
            + drag_force(CD0, CDa, alpha(vx, vz, pitch), alpha_i, rho, r, velocity(vx, y_velocity, vz))*-vz/velocity(vx, y_velocity, vz) ...
            - m*g)/m;
end
