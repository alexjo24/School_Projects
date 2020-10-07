function [ Ae_star ] = getAstar( lambda, Ae_min, Ae_max,l_0,q0e )
%Calculate A* for a given lambda. 

%Input: lambda
%       Ae_min, Minimum area
%       Ae_max, Maximum area
%       l_0 , Global length vector
%       q0e   , qo(e) (per element), computed in the main m-file.

%Output: A*

Ae_trial = sqrt(q0e/(lambda*l_0));


   if  Ae_trial < Ae_min;
       
  
          Ae_star = Ae_min;
          
   elseif Ae_trial > Ae_max 
       
       
       
         Ae_star = Ae_max;
   else
       
       
         Ae_star = Ae_trial;
   end



end

