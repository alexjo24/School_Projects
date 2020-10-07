function [ dphidlambda ] = dphidlambda_b(lambda,rho_max,rho_min,alpha,be,a_e,V_max,nelm,thickness)
%Compute dphidlambda, solving the dual function.
%
%Input: Lambda 
%       l_0 , Global length vector 
%       Vmax, Maximum volume
%       rho_min, Minimum density
%       rho_max, Maximum density
%       alpha, set parameter in the main m-file, connected to the
%       intervening variable in the SIMP OC method.
%       be   , b array, computed in the main m-file.
%       a_e  , area of each element
%       nelm , number of elements
%       thickness
%
%Output: dphidlambda

     h = 0;
    for z=1:nelm
        xe_star(z,:) = getx_e_star( lambda,rho_max,rho_min,alpha,be(z,:),a_e(z,:) );
       
        he = a_e(z,:)*xe_star(z,:);
        h = he + h;
    end 


%Constraint function
dphidlambda = h*thickness-V_max;

end

