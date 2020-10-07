function [ Ae_star ] = getAstar( lambda, Ae_min, Ae_max,l_0,q0e,Lj,alpha_j )
%Calculate A* for a given lambda. 

%Input: lambda
%       Ae_min, Minimum area
%       Ae_max, Maximum area
%       l_0 , Global length vector
%       q0e   , qo(e) (per element), computed in the main m-file.
%       Lj    , Moving asymptote, 
%       alpha_j, move limit

%Output: A*

Ae_trial =Lj + sqrt(q0e/(lambda*l_0));


   if  Ae_trial < alpha_j;
       
          Ae_star = alpha_j;
          
   elseif Ae_trial > Ae_max 
       
       
       
         Ae_star = Ae_max;
   else
       
       
         Ae_star = Ae_trial;
   end



end

