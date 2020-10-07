function [KinE, IntE] = plan3gEn( a_dot, M , D , E,nelm, ex,ey, th )
%Calculate the kinetic and internal energies of a three node element
% INPUT
% 
%  ex - global x-coordinates: [320x3]
%  ey - global y-coordinates  [320x3]
% 
% th = thickness
%
% a_dot - nodal velocity: [400x1]
%
% nelm - number of elements
% nelm - [320]
%
% M - Mass matrix
% M = [400x400]
%
% D = Constitutive matrix
% E = Element strain vector
%
% OUTPUT
% KinE - Kinetic energy
% IntE - Internal energy

IntE=0;

%Element loop 
for n=1:nelm
    
ec = [ex(n,:) ; ey(n,:)]; 

% Reference coordinates in the undeformed configuration        
x1_0 = ec(1,1);
x2_0 = ec(1,2);
x3_0 = ec(1,3);
y1_0 = ec(2,1);
y2_0 = ec(2,2);
y3_0 = ec(2,3);

%Compute Element Area.
A0 = (1/2)*abs(det([1 x1_0 y1_0;
                1 x2_0 y2_0;
                1 x3_0 y3_0]));

%Summation of internal energy regarding all elements at a certain load step
IntE = IntE+(1/2)*E(:,n)'*D*E(:,n)*th*A0;
end

%Compute kinetic energy
KinE = (1/2)*a_dot'*M*a_dot;


end

