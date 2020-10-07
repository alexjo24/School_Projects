function [ D ] = dmat1D( ep,eg )
%Calculating the 1D tangent stiffness

E = ep(1);
v = ep(2);


G = E/(2*(1+v));

% D = G*((1/sqrt(2*eg+1))+(2/((2*eg+1)^2)));
D = G*(1+2/((sqrt(2*eg+1))^3));


end

