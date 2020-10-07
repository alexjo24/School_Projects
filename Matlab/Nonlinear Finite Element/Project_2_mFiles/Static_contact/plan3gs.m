function [ ee,eff ] = plan3gs( ec,ed )
%Compute strains and deformation gradient in a triangular 3 node large
%deformation element.
%
% INPUT
% ec = [x1 x2;
%       y1 y2;
%       z1 z2]
%
% Nodal displacement
% ed = [u1 u2 u3 u4 u5 u6];
%
% OUTPUT
% ee - greens strain vector [3x1]
%
% eff - deformation gradien [4x1]

ec = ec';

x1_0 = ec(1,1);
x2_0 = ec(1,2);
x3_0 = ec(1,3);
y1_0 = ec(2,1);
y2_0 = ec(2,2);
y3_0 = ec(2,3);

%Compute element Area.
A0 = (1/2)*det([1 x1_0 y1_0;
                1 x2_0 y2_0;
                1 x3_0 y3_0]);
               
% Derivative of shapefunctions with respect of reference coordinates  
dN1_dx0 = (1/(2*A0))*(y2_0-y3_0);
dN2_dx0 = (1/(2*A0))*(y3_0-y1_0);
dN3_dx0 = (1/(2*A0))*(y1_0-y2_0);

dN1_dy0 = (1/(2*A0))*(x3_0-x2_0);
dN2_dy0 = (1/(2*A0))*(x1_0-x3_0);
dN3_dy0 = (1/(2*A0))*(x2_0-x1_0);


%Derivative of displacement with respect of reference coordinates,
% current load step
du1_dx0 = (ed(1)*dN1_dx0+ ed(3)*dN2_dx0+ ed(5)*dN3_dx0);
du2_dx0 = (ed(2)*dN1_dx0+ ed(4)*dN2_dx0+ ed(6)*dN3_dx0);

du1_dy0 = (ed(1)*dN1_dy0+ ed(3)*dN2_dy0+ ed(5)*dN3_dy0);
du2_dy0 = (ed(2)*dN1_dy0+ ed(4)*dN2_dy0+ ed(6)*dN3_dy0);

%Compute deformation gradient
F = eye(3) + [ du1_dx0 , du1_dy0 , 0;
               du2_dx0  , du2_dy0 , 0
                0       ,  0      , 0];


eff = double(F);

eff = [eff(1,1);
       eff(1,2);
       eff(2,1);
       eff(2,2)];

% Compute greens strain vector
C = F'*F;
ee = (1/2)*(C-eye(3));

ee = [ee(1,1);
      ee(2,2);
      2*ee(1,2)];
  

end

