m = 0;
x = linspace(0,1,20);
t = linspace(0,2,5);

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
% Extract the first solution component as u.  This is not necessary
% for a single equation, but makes a point about the form of the output.
u = sol(:,:,1);

% A surface plot is often a good way to study a solution.
figure;
surf(x,t,u);
title('Numerical solution computed with 20 mesh points.');
xlabel('Distance x');
ylabel('Time t');

figure;
surf(x,t,exp(-t)'*sin(pi*x));
title('True solution plotted with 20 mesh points.');
xlabel('Distance x');
ylabel('Time t');

% A solution profile can also be illuminating.
figure;
plot(x,u(end,:),'o',x,exp(-t(end))*sin(pi*x));
title('Solutions at t = 2.');
legend('Numerical, 20 mesh points','Analytical', 'Location', 'NorthEast');
xlabel('Distance x');
ylabel('u(x,2)');

% --------------------------------------------------------------------------

function [c,f,s] = pdex1pde(x,t,u,DuDx)
c = pi^2;
f = DuDx;
s = 0;

% --------------------------------------------------------------------------

function u0 = pdex1ic(x)
u0 = sin(pi*x);

% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = pi * exp(-t);
qr = 1;

function [CLo,CLa,al,ro,r,v]