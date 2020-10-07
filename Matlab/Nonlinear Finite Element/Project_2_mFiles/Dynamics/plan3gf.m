function [ ef ] = plan3gf( ec,t,ed,es )
%Compute internal element force vector in a triangular 3 node large
%deformation element in plane strain.
% 
% ec = [x1 x2 x3;
%      y1 y2 y3];
% 
%t = thickness
%
%ed = [a1 a2 a3 a4 a5 a6];
%
%ef = F_int' = [f1 f2 f3 f4 f5 f6];
%
%dx_dx_0 = diff(x,x0);
%
%
%

S = [es(1);
     es(2);
     es(3)];

x1_0 = ec(1,1);
x2_0 = ec(1,2);
x3_0 = ec(1,3);
y1_0 = ec(2,1);
y2_0 = ec(2,2);
y3_0 = ec(2,3);

%Element Area.
A0 = (1/2)*det([1 x1_0 y1_0;
                1 x2_0 y2_0;
                1 x3_0 y3_0]);
      
            
dN1_dx0 = (1/(2*A0))*(y2_0-y3_0);
dN2_dx0 = (1/(2*A0))*(y3_0-y1_0);
dN3_dx0 = (1/(2*A0))*(y1_0-y2_0);

dN1_dy0 = (1/(2*A0))*(x3_0-x2_0);
dN2_dy0 = (1/(2*A0))*(x1_0-x3_0);
dN3_dy0 = (1/(2*A0))*(x2_0-x1_0);            
                    

% u_e = [u1;
%        u2];     

du1_dx0 = (ed(1)*dN1_dx0+ ed(3)*dN2_dx0+ ed(5)*dN3_dx0);
du2_dx0 = (ed(2)*dN1_dx0+ ed(4)*dN2_dx0+ ed(6)*dN3_dx0);

du1_dy0 = (ed(1)*dN1_dy0+ ed(3)*dN2_dy0+ ed(5)*dN3_dy0);
du2_dy0 = (ed(2)*dN1_dy0+ ed(4)*dN2_dy0+ ed(6)*dN3_dy0);



B_0le = [dN1_dx0  ,  0      , dN2_dx0 , 0       , dN3_dx0,  0       ;
            0     , dN1_dy0  , 0       , dN2_dy0 , 0      , dN3_dy0  ;
            dN1_dy0,  dN1_dx0 , dN2_dy0 , dN2_dx0 , dN3_dy0 , dN3_dx0];
        
%eq (5.13)
Ae = [ du1_dx0   ,  0         , du2_dx0  ,     0   ;
          0      ,  du1_dy0   ,  0       ,  du2_dy0;
        du1_dy0  ,  du1_dx0   , du2_dy0  ,  du2_dx0];


H_0e = [dN1_dx0 , 0       , dN2_dx0  , 0       , dN3_dx0 ,    0      ;
        dN1_dy0 , 0       , dN2_dy0  , 0       , dN3_dy0 ,    0      ;
        0       , dN1_dx0 ,   0      , dN2_dx0 ,    0    , dN3_dx0   ;
        0       , dN1_dy0 ,   0      , dN2_dy0 ,    0    , dN3_dy0    ]; 
    


B_0 = B_0le + Ae*H_0e;

ef = B_0'*S*A0*t;



end

