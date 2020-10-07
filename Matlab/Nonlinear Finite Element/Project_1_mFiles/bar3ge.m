function [ Ke ] = bar3ge( ec,ep,ed,es )
%Calculating the stiffness matrix
%Ke = K0(x0)+Ku(x0,u_n-1)+K_sigma(epsilon n-1)
%
%ec = [x1 x2;
%      y1 y2;
%      z1 z2]
%
%ep = [E A0];
%ed = [u1 u2 u3 u4 u5 u6];

N = es;

EA = ep(1)*ep(2);

u = [ed(4)-ed(1);
    ed(5)-ed(2);
    ed(6)-ed(3)];

  x0 = [ec(1,2)-ec(1,1);
      ec(2,2)-ec(2,1);
      ec(3,2)-ec(3,1)];
   
  l_0 = sqrt(x0'*x0);
  
K_0 = ((EA / (l_0)^3))*[ x0*x0' -x0*x0'...
                        ; -x0*x0' x0*x0'];

K_u = ((EA / (l_0)^3))*[ (x0*u'+u*x0'+u*u') -(x0*u'+u*x0'+u*u') ....
                        ; -(x0*u'+u*x0'+u*u') (x0*u'+u*x0'+u*u')];

K_sigma = (N/l_0)*[eye(3)  -eye(3);...
                    -eye(3) eye(3)];
           
Ke = K_0 + K_u + K_sigma;
                    
end

