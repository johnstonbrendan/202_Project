function y = alpha(x_velocity,z_velocity, pitch)
y = atan(abs(z_velocity)/abs(x_velocity)) + pitch;
end