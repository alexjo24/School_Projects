function [ ef ] = bar3gf( ec,ed,es )
%Calculating the internal force vector for a bar element
%ec = [x1 x2;
%      y1 y2;
%      z1 z2]
%
%es = [N];
%
%ed = [u1 u2 u3 u4 u5 u6];


x0 = [ec(1,2)-ec(1,1);
      ec(2,2)-ec(2,1);
      ec(3,2)-ec(3,1)];

%Compute inital length of a bar element
l_0 = sqrt(x0'*x0); 

N = es;

x_A = [ec(1,1)+ed(1), ec(2,1)+ed(2), ec(3,1)+ed(3)];
x_B = [ec(1,2)+ed(4), ec(2,2)+ed(5), ec(3,2)+ed(6)];


x_tilde = [x_A , x_B];
ef = (N/l_0)*[eye(3) -eye(3); -eye(3) eye(3)]*x_tilde';

ef = ef';



end

