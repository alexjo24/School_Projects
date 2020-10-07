function [ phi_lambda ] = phi_lambda( lambda,rho_max,rho_min,alpha,be,a_e,V_max,nelm )
%Compute dphidlambda

     phi_e = 0;
    for z=1:nelm
        xe_star(z,:) = getx_e_star( lambda,rho_max,rho_min,alpha,be(z,:),a_e(z,:) );
        
       phi_e = phi_e + be(z,:)*(xe_star(z,:)^(-alpha))+lambda*a_e(z,:)*xe_star(z,:);
    end 


  phi_lambda= phi_e -lambda*V_max;

end

