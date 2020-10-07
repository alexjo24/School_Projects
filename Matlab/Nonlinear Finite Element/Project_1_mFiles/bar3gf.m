function [ ef ] = bar3gf( ec,ed,es )
%Calculating the internal force vector
%ec = [x1 x2;
%      y1 y2;
%      z1 z2]
%es = [N];
%
%ed = [a1 a2 a3 a4 a5 a6];


x0 = [ec(1,2)-ec(1,1);
      ec(2,2)-ec(2,1);
      ec(3,2)-ec(3,1)];
    
 l_0 = sqrt(x0'*x0); 


N = es;

x_A = [ec(1,1)+ed(1), ec(2,1)+ed(2), ec(3,1)+ed(3)];
x_B = [ec(1,2)+ed(4), ec(2,2)+ed(5), ec(3,2)+ed(6)];


%page 34, q_tilde
x_tilde = [x_A , x_B];
ef = (N/l_0)*[eye(3) -eye(3); -eye(3) eye(3)]*x_tilde';

ef = ef';



end

