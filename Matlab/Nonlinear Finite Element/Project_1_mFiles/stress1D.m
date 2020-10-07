function [S] = stress1D(E,v,ec,ed)
%Calculating the secon Piola-Kirchhoff stress
%     x0_1=ec(1,1);
%     x0_2=ec(1,2);
%     y0_1=ec(2,1);
%     y0_2=ec(2,2);
%     z0_1=ec(3,1);
%     z0_2=ec(3,2);
    
    x_0=[ec(1,2)-ec(1,1);
         ec(2,2)-ec(2,1);
         ec(3,2)-ec(3,1)];
    
      x=[(ec(1,2)+ed(4))-(ec(1,1)+ed(1));
         (ec(2,2)+ed(5))-(ec(2,1)+ed(2));
         (ec(3,2)+ed(6))-(ec(3,1)+ed(3))];
    
l_0=sqrt(x_0'*x_0);
l=sqrt(x'*x);

LAMBDA=l/l_0;
Ges=E/(2*(1+v));
    
S=Ges*(LAMBDA-1/LAMBDA^2);
% S=S*10^6;

end

