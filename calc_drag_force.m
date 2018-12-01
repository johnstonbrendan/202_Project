function y = calc_drag_force(CD0, CDa, a, alpha_i, rho, r, v)
y = 1/2*(CD0+CDa*(a-alpha_i)^2)*rho*pi*(r^2)*(v^2); 
end