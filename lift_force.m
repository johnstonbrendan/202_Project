function y = lift_force(CL0,CLa,alpha,rho,r,v)
y = 1/2*(CL0+CLa*alpha)*rho*pi*power(r,2)*power(v,2);
end