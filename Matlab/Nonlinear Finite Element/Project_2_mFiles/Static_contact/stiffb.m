function [D] = stiffb(ec,ee,k,rc) 
%Calculate material stiffness for a bar
%Penalty method is implemented
%
% INPUT
% ec = [x1 x2;
%       y1 y2;
%       z1 z2]
%
% ee = strain of a bar element
%
%k = stiffness of the bar element
%
%rc = radius of cylinder
%
% OUTPUT
% D - tangent stiffness tensor

% Reference coordinates in the undeformed configuration
x0 = [ec(1,2)-ec(1,1);
      ec(2,2)-ec(2,1);
      ec(3,2)-ec(3,1)];
    

% Initial length of the bar     
l_0 = sqrt(x0'*x0); 

        
 %Stretch as a function of greens strain        
 Lambda = sqrt(2*ee+1);

 Lambda_c = rc/l_0;  

    %Constitutive law, Penalty method regarding bar elements 
    if Lambda < Lambda_c
        
       %Tangent stiffness tensor
       D = k*Lambda_c*(Lambda^(-3));
        
    else
        %Tangent stiffness tensor
        D = 0;
        
    end

end

