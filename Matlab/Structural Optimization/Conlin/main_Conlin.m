clc
clear all
close all
format long
tic

load geomSO.mat

%Material parameters
E = 1;      %Young's modulus
Ae_max = (((20*10^-3)^2)*pi)/4;          %Maximum cross sectional area
Ae_min = (((20*10^-8)^2)*pi)/4;          %Minimum cross sectional area

V_max = (2*10^-6)/2; % Maximum volume allowed


 %Number of elements
 nelm = max(edof(:,1)); 
 
 %Initiate 
 K = sparse(ndof,ndof);           %Introduce stiffness matrix, Ke [8x8]
 A_global = ones(nelm,1)*1*10^-8; %Introduce element Area vector, with all Area being

%Initate Global cell for normalized stiffness matrix
 Cell_K0e = cell(nelm,1);
 
    for n=1:nelm
        
        
     %Initate constant norm element stiffness matrix
     Cell_K0e{n} = bar2e(ex(n,:),ey(n,:),[E 1]);     
     
     %Obtain a global length vector for each element.
      ec =[ex(n,:),ey(n,:)];
      
     x_0=[ec(1,2) - ec(1,1);
          ec(1,4) - ec(1,3)];

      l_0 = sqrt(x_0'*x_0);
        
      
     %Global length vector
     l_0_global(n,:) = l_0;
 
    end

    
 V = A_global'*l_0_global;   
    
 %Initate Tolerance
 Tol = 10^-9;
 
 %Initate quantaties
 iter = 0;
 Nres = 0;
 q0 = 0;
 
 %Input min and max values for fzero function.
 lambda_min=10^10;
 lambda_max=10^15;
 
  
%Optimization loop
while Nres > Tol || iter == 0
    
    iter = iter + 1;   
    disp(['Iteration ' ,num2str(iter),'----------------------------'])
    
    
    %Reset global stiffness matrix
     K=sparse(ndof,ndof);
      
    %Element loop to assemble global stiffness matrix
      for i=1:nelm

          ep = [E A_global(i,:)];

          Ke = bar2e(ex(i,:),ey(i,:),ep);

         %Assemble element matrices to global matrices:
             indx = edof(i,2:end);
             K(indx,indx) = K(indx,indx)+Ke;
      end

      %Obtain global displacements
          a = solveq(K,f,bc);
          ed = extract(edof,a);
         
    for j = 1:nelm
            
      %Compute q0j according to eq 5.37, applied for CONLIN approx
      q0e(j,:) = (ed(j,:)*Cell_K0e{j}*ed(j,:)'*(A_global(j,:)^2));

    end

        
        %Solve the dual function
        lambdastar = fzero(@(lambda) dphidlambda...
        ( lambda,l_0_global,V_max,Ae_min, Ae_max,nelm,q0e )...
        ,[lambda_min lambda_max]);
    

for z=1:nelm

       %Find new rho value per element based on new lambda from the dual
       %function
       Ae_star(z,:)= getAstar( lambdastar , Ae_min, Ae_max,...
                                l_0_global(z,:),q0e(z,:) );
       
       %Compute element forces
       stress(z,:) = bar2s(ex(z,:),ey(z,:),[E Ae_star(z,:)],ed(z,:))/Ae_star(z,:);
  
end       
     
     K=sparse(ndof,ndof);
     %Element loop to assemble global stiffness matrix
      for i=1:nelm

          ep = [E A_global(i,:)];

          Ke = bar2e(ex(i,:),ey(i,:),ep);


         %Assemble element matrices to global matrices:
             indx = edof(i,2:end);
             K(indx,indx) = K(indx,indx)+Ke;
      end

      %Obtain displacement
          a = solveq(K,f,bc);
          ed = extract(edof,a);
          
 for z=1:nelm

     %Compute element stress
     stress(z,:) = bar2s(ex(z,:),ey(z,:),[E Ae_star(z,:)],ed(z,:))/Ae_star(z,:);
     
 end
          
    %Compute the residual between current and previous iteration
    %regarding the updated area.
     Nres = norm(Ae_star-A_global);
     A_global = Ae_star;
     
     %Plot vectors
     Ae_star_plot(:,iter) = Ae_star;
     Iter_plot(:,iter) = iter;
     g0(:,iter) = f'*a;
     stress_global(:,iter) = stress;
    
     g1(:,iter) = A_global'*l_0_global-V_max;

     V = A_global'*l_0_global;
     
     fprintf('Iter %i Res %e\n',iter,Nres)
     
     end

%%

figure
fac = 4*10^6;      %Magnification of Cross Section area.
magnfac = 1*10^-10; %Magnification of Deformed structure.
myeldisp2(ex,ey,ed,[1 2 0],magnfac,A_global,fac);
grid on
xlabel('x,[m]')
ylabel('y,[m]')
title(['Structural optimization with CONLIN approximation'],'FontSize', 14)

figure
plot(Iter_plot,g0,'LineWidth',2)
hold on

plot(Iter_plot,g1,'LineWidth',2)
grid on
xlabel('Iterationer','FontSize', 14)
ylabel('g_i','FontSize', 14)
legend('g0- objective function','g1- constraint function')
title('g_i functions with CONLIN approximation','FontSize', 14)
set(gca,'FontSize',14)


%Check Area fulfills Box constraint

A_max_check = max(A_global);
A_min_check = min(A_global);

    if A_max_check <= Ae_max
        disp('Valid Box constraint regarding A_max')
    else
        disp('Box constraint not fulfilled')
    end

    
     if A_min_check >= Ae_min
        disp('Valid Box constraint regarding A_min')
        else
        disp('Box constraint not fulfilled')
     end
 toc
  %Check procentage of elements in A_global which have the values Ae_min
  %and Ae_max.
  procent_min=length(find(A_global == Ae_min))/nelm 
  procent_max=length(find(A_global == Ae_max))/nelm 
  
  