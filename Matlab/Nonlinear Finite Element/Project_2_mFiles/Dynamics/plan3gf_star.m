function [ ef ] = plan3gf_star( ec,th,ed_n,ed,stress ,stress_n )
%Compute internal element force vector in a triangular 3 node large
%deformation element in plane strain. Mean values are computed over last
%equilibrium load step and current load step.
%
% INPUT
% ec = [x1 x2 x3;
%      y1 y2 y3];
% 
% th = thickness
%
% Nodal displacement
% ed = [u1 u2 u3 u4 u5 u6];
%
% Stress vector
% stress = [S11
%           S22
%           S12];
%
% OUTPUT
% ef=Fint'=[f1, f2,.., f6]

%Stress vector
S = [stress(1);
     stress(2);
     stress(3)];

 % Reference coordinates in the undeformed configuration
x1_0 = ec(1,1);
x2_0 = ec(1,2);
x3_0 = ec(1,3);
y1_0 = ec(2,1);
y2_0 = ec(2,2);
y3_0 = ec(2,3);

%Compute Element Area.
A0 = (1/2)*det([1 x1_0 y1_0;
                1 x2_0 y2_0;
                1 x3_0 y3_0]);
      
%Derivative of shapefunctions with respect of reference coordinates            
dN1_dx0 = (1/(2*A0))*(y2_0-y3_0);
dN2_dx0 = (1/(2*A0))*(y3_0-y1_0);
dN3_dx0 = (1/(2*A0))*(y1_0-y2_0);

dN1_dy0 = (1/(2*A0))*(x3_0-x2_0);
dN2_dy0 = (1/(2*A0))*(x1_0-x3_0);
dN3_dy0 = (1/(2*A0))*(x2_0-x1_0);            
                       

%Derivative of displacement with respect of reference coordinates
% last equilibrum step
du1_dx0_n = (ed_n(1)*dN1_dx0+ ed_n(3)*dN2_dx0+ ed_n(5)*dN3_dx0);
du2_dx0_n = (ed_n(2)*dN1_dx0+ ed_n(4)*dN2_dx0+ ed_n(6)*dN3_dx0);

du1_dy0_n = (ed_n(1)*dN1_dy0+ ed_n(3)*dN2_dy0+ ed_n(5)*dN3_dy0);
du2_dy0_n = (ed_n(2)*dN1_dy0+ ed_n(4)*dN2_dy0+ ed_n(6)*dN3_dy0);


%Derivative of displacement with respect of reference coordinates,
% current load step
du1_dx0_n1 = (ed(1)*dN1_dx0+ ed(3)*dN2_dx0+ ed(5)*dN3_dx0);
du2_dx0_n1 = (ed(2)*dN1_dx0+ ed(4)*dN2_dx0+ ed(6)*dN3_dx0);

du1_dy0_n1 = (ed(1)*dN1_dy0+ ed(3)*dN2_dy0+ ed(5)*dN3_dy0);
du2_dy0_n1 = (ed(2)*dN1_dy0+ ed(4)*dN2_dy0+ ed(6)*dN3_dy0);



B_0le = [dN1_dx0  ,  0      , dN2_dx0 , 0       , dN3_dx0,  0       ;
            0     , dN1_dy0  , 0       , dN2_dy0 , 0      , dN3_dy0  ;
            dN1_dy0,  dN1_dx0 , dN2_dy0 , dN2_dx0 , dN3_dy0 , dN3_dx0];
        

Ae_n = [ du1_dx0_n   ,  0         , du2_dx0_n  ,     0   ;
          0      ,  du1_dy0_n   ,  0       ,  du2_dy0_n;
        du1_dy0_n  ,  du1_dx0_n   , du2_dy0_n  ,  du2_dx0_n];
    
    
    
Ae_n1 = [ du1_dx0_n1   ,  0         , du2_dx0_n1  ,     0   ;
          0      ,  du1_dy0_n1   ,  0       ,  du2_dy0_n1;
        du1_dy0_n1  ,  du1_dx0_n1   , du2_dy0_n1  ,  du2_dx0_n1];
    
Ae = (Ae_n1 + Ae_n)*(1/2);

H_0e = [dN1_dx0 , 0       , dN2_dx0  , 0       , dN3_dx0 ,    0      ;
        dN1_dy0 , 0       , dN2_dy0  , 0       , dN3_dy0 ,    0      ;
        0       , dN1_dx0 ,   0      , dN2_dx0 ,    0    , dN3_dx0   ;
        0       , dN1_dy0 ,   0      , dN2_dy0 ,    0    , dN3_dy0    ]; 
    


B_0 = B_0le + Ae*H_0e;

% Mean stress matrix
S_star = (1/2)*(stress + stress_n);


%Compute element internal force vector
ef = B_0'*S_star*A0*th;



end

