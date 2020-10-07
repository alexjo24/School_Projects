function [ dphidlambda ] = dphidlambda_b(lambda,rho_max,rho_min,alpha,be,a_e,V_max,nelm,th)
%Compute dphidlambda

     h = 0;
    for z=1:nelm
        xe_star(z,:) = getx_e_star( lambda,rho_max,rho_min,alpha,be(z,:),a_e(z,:) );
       
        he = a_e(z,:)*xe_star(z,:);
        h = he + h;
    end 



dphidlambda = h*th-V_max;

end

