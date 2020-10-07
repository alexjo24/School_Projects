clc
clear all
close all

%% Manually set variables
%Number of grid points manually entered by user
nbr_of_gridpoints = 18;
%Number of iterations manually entered by user
nbr_iter = 200;

%%
%10 node point mesh
Grid = 0:nbr_of_gridpoints; 

%Length of nozzle is 3L
L = 0.05;

x = (3*L/(length(Grid)-1))*Grid; 

x_prim = x./L;

A_prim =  1+2.2*((x_prim-1.5).^2);



%Initial conditions
rho_prim = 1-0.3146*x_prim;
T_prim = 1-0.2314*x_prim;
V_prim = (0.1+1.09*x_prim).*sqrt(T_prim);


%Local speed of sound
gamma = 1.4; %Standard air conditions
R = 287.058; %Universal gas constant for dry air 
% a = sqrt(gamma*R*T_prim);
a = sqrt(T_prim);

%Inviscid stability condition

%Conservative approach, set Mach = 3 for stability condition regarding
%Maximum possible velocities in the system for all iterations.


% M__stability = 2; %Mach number 
% T_stability = 500;  %Temperature
% a_stability = sqrt(gamma*R*T_stability);
% V_stability = M__stability*a_stability;
% delta_x = x_prim(2)-x_prim(1);
% delta_t = ((V_stability/delta_x)+a_stability*sqrt((1/((delta_x)^2))))^-1;
% dtdx = delta_t/delta_x;

%Courant-Friedrichs-Lewy condition
a_prim = sqrt(T_prim);
delta_x = x_prim(2)-x_prim(1);
C = 0.8;
delta_t = min(C*(delta_x./(a+V_prim)));

dtdx = delta_t./delta_x;


%"The residual is the difference between the previous result and the 
%current result. As these errors are decreasing the equation results 
%are reaching values that are changing less and less. This is what is 
%known as convergence. That is the solutions are converging."
%Reference: https://www.cfd-online.com/Forums/fluent/28009-residuals.html

%I decide to use a for loop for manual control of number of iterations,
%while loop may also be used.

%% Boundary conditions version 1

%Subsonic inlet
% M < 0.8
% 1 Mach = Velocity / local speed of sound

%In Isentropic flow properties table
%Mach = 0.6 provides:
%                    T0 / T = 1.072
%                    p0 / p = 1.276
%                    rho0 / rho = 1.19
%                    A / A_star = 1.18


% M_inlet = 0.6; %Inlet Mach number 
% T_inlet = 400;  %Inflow temp 400 K
% a_inlet = sqrt(gamma*R*T_inlet);
% 
% T0_inlet = T_inlet*1.072; 
% a0_inlet = sqrt(gamma*R*T0_inlet);
% 
% V_inlet = M_inlet*a_inlet;
% V_prim_inlet = V_inlet/a0_inlet;
% T_prim_inlet = T_inlet/T0_inlet;
% rho_prim_inlet = 1/1.19;
% 
% 
% %Set inlet boundary conditions
%     V_prim(1) = V_prim_inlet;
%     T_prim(1) = T_prim_inlet;
%     rho_prim(1) = rho_prim_inlet;


%% Boundary conditions version 2

        %Extrapolation at boundaries if needed
        %Inlet: p1 = 2*p2 - p3
        %Outlet: pN = 2*pNmin1 - pNmin2   
        
        %Subsonic inflow
        %Extrapolate 1 boundary condition from interior

        V_prim_inlet_extrapo =  2*V_prim(2)-V_prim(3);
        
        
        %Supersonic outflow
        %Extrapolate all boundary conditions from interior
        
        
        rho_prim_outlet_extrapo =  2*rho_prim(length(rho_prim)-1)-...
                                  rho_prim(length(rho_prim)-2);
                              
        V_prim_outlet_extrapo =  2*V_prim(length(V_prim)-1)-...
                                  V_prim(length(V_prim)-2);
                              
        T_prim_outlet_extrapo =  2*T_prim(length(T_prim)-1)-...
                                  T_prim(length(T_prim)-2);

rho_prim_inlet = 1;
T_prim_inlet = 1;

rho_prim(1) = rho_prim_inlet;
rho_prim(length(rho_prim)) = rho_prim_outlet_extrapo;

T_prim(1) = T_prim_inlet;
T_prim(length(T_prim)) = T_prim_outlet_extrapo;

V_prim(1) = V_prim_inlet_extrapo;
V_prim(length(V_prim)) = V_prim_outlet_extrapo;




%% Start of iteration loop
for iter = 1:nbr_iter

 
    
for j = 2:(length(Grid)-1)

%A' values does not alternate between timestep so it is possible to compute
% outside the for loop
dlnA_dx_ver1 = log(A_prim(j+1))-log(A_prim(j));    
dlnA_dx_ver2 = log(A_prim(j))-log(A_prim(j-1));    


%% Predictor
%Continuity

rho_nplus1_j_s(j-1) = rho_prim(j)-V_prim(j)*dtdx*...
                    (rho_prim(j+1)-rho_prim(j))-rho_prim(j)*dtdx*...
                    (V_prim(j+1)-V_prim(j))-rho_prim(j)*V_prim(j)*...
                    dtdx*dlnA_dx_ver1; 
 
                
                
rho_nplus1_jmin1_s(j-1) = rho_prim(j-1)-V_prim(j-1)*dtdx*...
                    (rho_prim(j)-rho_prim(j-1))-rho_prim(j-1)*dtdx*...
                    (V_prim(j)-V_prim(j-1))-rho_prim(j-1)*V_prim(j-1)*...
                    dtdx*dlnA_dx_ver2;                 

%Momentum

V_nplus1_j_s(j-1) = V_prim(j)- V_prim(j)*dtdx*(V_prim(j+1)-V_prim(j))-...
                  (1/gamma)*dtdx*((T_prim(j+1)-T_prim(j))+...
                  (T_prim(j)/rho_prim(j))*(rho_prim(j+1)-rho_prim(j)));
              
V_nplus1_jmin1_s(j-1) = V_prim(j-1)- V_prim(j-1)*dtdx*(V_prim(j)-...
                V_prim(j-1))-(1/gamma)*dtdx*((T_prim(j)-T_prim(j-1))+...
                  (T_prim(j-1)/rho_prim(j-1))*(rho_prim(j)-rho_prim(j-1)));
              
%Energy              
T_nplus1_j_s(j-1) = T_prim(j)-V_prim(j)*dtdx*(T_prim(j+1)-T_prim(j))-...
                  (gamma-1)*T_prim(j)*dtdx*((V_prim(j+1)-V_prim(j))+...
                  V_prim(j)*dlnA_dx_ver1);              
              
T_nplus1_jmin1_s(j-1) = T_prim(j-1)-V_prim(j-1)*dtdx*(T_prim(j)-...
                  T_prim(j-1))-(gamma-1)*T_prim(j-1)*dtdx*...
                  ((V_prim(j)-V_prim(j-1))+...
                  V_prim(j-1)*dlnA_dx_ver2);
              
%% Corrector
%Continuity

rho_nplus1_j(j-1) = 0.5*((rho_prim(j)+rho_nplus1_j_s(j-1)) - ...
                    dtdx*V_nplus1_j_s(j-1)*(rho_nplus1_j_s(j-1)-...
                    rho_nplus1_jmin1_s(j-1))-dtdx*rho_nplus1_j_s(j-1)*...
                    (V_nplus1_j_s(j-1)-V_nplus1_jmin1_s(j-1))-dtdx*...
                    rho_nplus1_j_s(j-1)*V_nplus1_j_s(j-1)*dlnA_dx_ver2);
                
                
V_nplus1_j(j-1) = 0.5*((V_prim(j)+V_nplus1_j_s(j-1))-dtdx*V_nplus1_j_s(j-1)*...
                (V_nplus1_j_s(j-1)-V_nplus1_jmin1_s(j-1))-...
               dtdx*(1/gamma)*((T_nplus1_j_s(j-1)-T_nplus1_jmin1_s(j-1))...
                +(T_nplus1_j_s(j-1)/rho_nplus1_j_s(j-1))*(rho_nplus1_j_s(j-1)...
                -rho_nplus1_jmin1_s(j-1))));
            
            
T_nplus1_j(j-1) = 0.5*((T_prim(j)+T_nplus1_j_s(j-1))-dtdx*V_nplus1_j_s(j-1)*...
                (T_nplus1_j_s(j-1)-T_nplus1_jmin1_s(j-1))-...
                (gamma-1)*T_nplus1_j_s(j-1)*dtdx*((V_nplus1_j_s(j-1)-...
                V_nplus1_jmin1_s(j-1))+V_nplus1_j_s(j-1)*dlnA_dx_ver2));


end



%Check residuals

rho_res = norm(rho_prim(2:(length(Grid)-1))-rho_nplus1_j);
V_res = norm(V_prim(2:(length(Grid)-1))-V_nplus1_j);
T_res = norm(T_prim(2:(length(Grid)-1))-T_nplus1_j);
disp('----------------------------------------------')
disp([num2str(iter),' Rho Residual ',num2str(rho_res)])
disp([num2str(iter),' Velocity Residual ',num2str(V_res)])
disp([num2str(iter),' Temperature Residual ',num2str(T_res)])


rho_prim = rho_nplus1_j;
T_prim = T_nplus1_j;
V_prim = V_nplus1_j;

        %Extrapolation at boundaries if needed
        %Inlet: p1 = 2*p2 - p3
        %Outlet: pN = 2*pNmin1 - pNmin2   
        
        %Subsonic inflow
        %Extrapolate 1 boundary condition from interior
        
        V_prim_inlet_extrapo =  2*V_prim(1)-V_prim(2);
        
        
        %Supersonic outflow
        %Extrapolate all boundary conditions from interior
        
        
        rho_prim_outlet_extrapo =  2*rho_prim(length(rho_prim))-...
                                  rho_prim(length(rho_prim)-1);
                              
        V_prim_outlet_extrapo =  2*V_prim(length(V_prim))-...
                                  V_prim(length(V_prim)-1);
                              
        T_prim_outlet_extrapo =  2*T_prim(length(T_prim))-...
                                  T_prim(length(T_prim)-1);


rho_prim = [rho_prim_inlet rho_prim rho_prim_outlet_extrapo];
T_prim = [T_prim_inlet T_prim T_prim_outlet_extrapo];
V_prim = [V_prim_inlet_extrapo V_prim V_prim_outlet_extrapo];

%Note the rest of the zero values will obtain a value in the next 
%iteration from extrapolation

%Arrays for plot
if iter == 1
    iter_arr = iter;
    rho_res_arr = rho_res; 
    V_res_arr = V_res;
    T_res_arr = T_res;
else
    iter_arr = [iter_arr iter];
    rho_res_arr = [rho_res_arr rho_res];  
    V_res_arr = [V_res_arr V_res];
    T_res_arr = [T_res_arr T_res];
end


clf
%Plot of residuals
plot(iter_arr,rho_res_arr);
grid on
hold on
plot(iter_arr,V_res_arr)
hold on
plot(iter_arr,T_res_arr)
% title(['Load step:', num2str(n)],'FontSize', 20)
title('Residual plot','FontSize', 20)
xlabel('Iterations','FontSize', 14)
ylabel('Residuals','FontSize', 14)
set(gca,'FontSize',14)
legend('Rho','Velocity','Temperature')
drawnow



end


%% Plot

%Nondimensional Mach number
M_plot = V_prim.*sqrt(abs(T_prim));

figure
plot(Grid,V_prim)
grid on
title('Velocity distribution')
xlabel('Grid points')
ylabel('Non-dimensional velocity')

figure
bar(V_prim)
title('V prim')

figure
bar(A_prim)
title('A prim')

figure
bar(T_prim)
title('T prim')

figure
bar(rho_prim)
title('rho prim')

figure
bar(M_plot)
title('M plot')