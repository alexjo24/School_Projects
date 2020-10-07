function [ xe_star ] = getx_e_star( lambda,rho_max,rho_min,alpha,be,a_e )
%Calculate A* for a given lambda. 


xe_trial = ((alpha*be)/(lambda*a_e))^(1/(1+alpha));


   if  xe_trial < rho_min;
       
  
          xe_star = rho_min;
          
   elseif xe_trial > rho_max 
       
       
       
         xe_star = rho_max;
   else
       
       
         xe_star = xe_trial;
   end



end

