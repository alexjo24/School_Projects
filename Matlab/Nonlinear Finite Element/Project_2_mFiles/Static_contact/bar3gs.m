function [  ee ] = bar3gs( ec,ed )
%[ es, ee ] = bar3gs( ec,ep,ed )
%Calculating Green's strain and the corresponding normal force
%for a bar element
%ec = [x1 x2;
%      y1 y2;
%      z1 z2];
%
%Nodal displacement
%ed = [u1 u2 u3 u4 u5 u6];

x0 = [ec(1,2)-ec(1,1);
      ec(2,2)-ec(2,1);
      ec(3,2)-ec(3,1)];



x_A = [ec(1,1)+ed(1), ec(2,1)+ed(2), ec(3,1)+ed(3)];
x_B = [ec(1,2)+ed(4), ec(2,2)+ed(5), ec(3,2)+ed(6)];


x_A_0 = [ec(1,1), ec(2,1), ec(3,1)];
x_B_0 = [ec(1,2), ec(2,2), ec(3,2)];


%Node A and node B for one bar element.
x_tilde = [x_A , x_B];

x_0_tilde = [x_A_0, x_B_0];


%Compute inital length of a bar element
l_0 = sqrt(x0'*x0);

       
u_tilde = [ed(1) ed(2) ed(3) ed(4) ed(5) ed(6)];


% Green's strain
ee = (1/((l_0)^2))*(1/2)*(x_0_tilde+x_tilde)*[eye(3) -eye(3); -eye(3) eye(3)]*u_tilde';


end

