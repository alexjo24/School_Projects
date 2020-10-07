function [ xe_star ] = getx_e_star( lambda,rho_max,rho_min,alpha,be,a_e )
%Calculate x* (= rho*) for a given lambda. 

%Input: lambda
%       rho_min, Minimum density
%       rho_max, Maximum density
%       l_0 , Global length vector
%       be   , b array, computed in the main m-file.
%       alpha, set parameter in the main m-file, connected to the
%       intervening variable in the SIMP OC method.
%       a_e  , area of each element

%Output: x*


xe_trial = ((alpha*be)/(lambda*a_e))^(1/(1+alpha));


   if  xe_trial < rho_min;
       
  
          xe_star = rho_min;
          
   elseif xe_trial > rho_max 
       
       
       
         xe_star = rho_max;
   else
       
       
         xe_star = xe_trial;
   end



end

