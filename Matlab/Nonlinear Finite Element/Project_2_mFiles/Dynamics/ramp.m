function [ P, delta_t ] = ramp( loadsteps )
%Ramp function

delta_t = linspace(0,(10^-3),loadsteps)';

P=(delta_t.^2);


end

