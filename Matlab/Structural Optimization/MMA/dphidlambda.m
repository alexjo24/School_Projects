function [ dphidlambda ] = dphidlambda( lambda,l_0_global,Vmax,Ae_min, Ae_max,nelm,q0e,Lj,alpha_j )
%Compute dphidlambda, solving the dual function.
%
%Input: Lambda 
%       l_0 , Global length vector 
%       Vmax, Maximum volume
%       Ae_min, Minimum area
%       Ae_max, Maximum area
%       q0e   , qo array, computed in the main m-file.
%       Lj    , Moving asymptote, 
%       alpha_j, move limit

%Output: dphidlambda

     h = 0;
    for z=1:nelm
        Ae_star(z,:) = getAstar( lambda , Ae_min, Ae_max,l_0_global(z,:),q0e(z,:),Lj(z,:),alpha_j(z,:) );
       
        he = l_0_global(z,:)*Ae_star(z,:);
        h = he + h;
    end 


%Constraint function
dphidlambda = h-Vmax;

end

