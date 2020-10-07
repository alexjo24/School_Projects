function [ Ke ] = plan3ge( ec,t,D,ed,es )
%Compute the element stiffness matrix for a triangular 3 node large
%deformation element in plane strain
%
%ec = [x1 x2 x3;
%      y1 y2 y3];
%
%ed = [u1 u2 u3 u4 u5 u6];    
%
%D = [3x3 matrix];
%
%es = [S11;
%      S22;
%      S12];
%

ec = ec';

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

S_matrix = [es(1) , es(3);
            es(3) , es(2)];
        
R = [ S_matrix , zeros(2);
      zeros(2) , S_matrix];

Ke = B_0'*D*B_0*(A0*t)+H_0e'*R*H_0e*(A0*t);


                    
end

