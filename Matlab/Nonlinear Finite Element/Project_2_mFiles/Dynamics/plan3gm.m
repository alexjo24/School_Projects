function [ Me ] = plan3gm( ec,th,rho )
%Compute the mass matrix for a 3 node element
% INPUT
%ec = [x1 x2 x3;
%      y1 y2 y3];
%
% rho - density
%
% th -thickness
%
% OUTPUT
% Me - Element mass matrix
%



% Reference coordinates in the undeformed configurationx1_0 = ec(1,1);
x1_0 = ec(1,1);
x2_0 = ec(1,2);
x3_0 = ec(1,3);
y1_0 = ec(2,1);
y2_0 = ec(2,2);
y3_0 = ec(2,3);


%Dunavant rule used, 
rule = 2;

%Three integration points used.
[ node_xy, w ] = dunavant_rule( rule);
w = 0.5*w;
node_xy = [node_xy(1,2) node_xy(1,1) node_xy(1,3);
           node_xy(2,2) node_xy(2,1) node_xy(2,3)];



%Compute Jacobian
J = [-x1_0+x2_0  -x1_0+x3_0
     -y1_0+y2_0   -y1_0+y3_0];


%Initiate mass matrix
Me = 0;

for n=1:3
    
    %Shape functions
    N1 = 1-node_xy(1,n)-node_xy(2,n);
    N2 = node_xy(1,n);
    N3 = node_xy(2,n);
    
    N = [N1 0 N2 0 N3 0 
         0  N1 0 N2 0 N3];
   
     
  %Assemble global mass matrix
  Me = Me + rho*th*det(J)*w(n).*N'*N;
end


end

