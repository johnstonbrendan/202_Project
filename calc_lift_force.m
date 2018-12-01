function y = calc_lift_force(CL0,CLa,a,rho,r,v)
y = 1/2*(CL0+CLa*a)*rho*pi*power(r,2)*power(v,2);
end