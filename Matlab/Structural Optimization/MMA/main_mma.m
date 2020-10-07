clc
clear all
close all
format long

tic

load geomSO.mat

%Material parameters
E = 1;      %Young's modulus
Ae_max = (((20*10^-3)^2)*pi)/4;     %Maximum cross sectional area
Ae_min = (((20*10^-8)^2)*pi)/4;    %Minimum cross sectional area
V_max = 2*10^-6;                    %Maximum allowed volume


%Number of elements
nelm = max(edof(:,1));

%Initiate
K = sparse(ndof,ndof);
A_global = ones(nelm,1)*10^-8; %Introduce element Area vector, with all Area being

%Initate Global cell for normalized stiffness matrix
Cell_K0e = cell(nelm,1);


%Initiate MMA parameters
s_slower = 0.6;
s_faster = 1.1;
s_init = 0.1;
A_last = 0;
A_lastlast = 0;
my = 0.5;


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


lambda_min=10^10;
lambda_max=10^15;


%Optimization loop
while Nres > Tol || iter == 0
    
    iter = iter + 1;
    disp(['Iteration ' ,num2str(iter),'----------------------------'])
    
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
    
    for j=1:nelm
        
        
        %Modifying the asymptotes during the iterations
        %according to page 68 in the course literature.
        if iter < 3
            Lj(j,:) = A_global(j,:)-s_init*(Ae_max-Ae_min);
        else
            
            sign1 = A_global(j,:) - A_last(j,:);
            sign2 = A_last(j,:) - A_lastlast(j,:);
            
            if sign1*sign2 > 0
                Lj(j,:) = A_global(j,:) - s_faster*(A_last(j,:)-Lj(j,:));
                
            else
                Lj(j,:) = A_global(j,:)-s_slower*(A_last(j,:)-Lj(j,:));
                
            end
        end
        
        
        %If Lj is negative, give it a value of zero.
        if sign(Lj(j,:)) ~= 1
            Lj(j,:) = 0;
            %disp('Negative Lj')
        end
        
        %Update alpha_j per element
        alpha_j(j,:) = max([Ae_min, Lj(j,:)+my*(A_global(j,:)-Lj(j,:))]);
        
        %Compute q0j according to eq 5.37
        q0e(j,:) =((A_global(j,:)-Lj(j,:))^2)*ed(j,:)*Cell_K0e{j}*ed(j,:)';
        
        
        
    end
    
    %Solve the dual function
    lambdastar = fzero(@(lambda) dphidlambda...
        ( lambda,l_0_global,V_max,Ae_min, Ae_max,nelm,q0e, Lj, alpha_j )...
        ,[lambda_min lambda_max]);
    
    for z=1:nelm
        
        %Find new rho value per element based on new lambda from the dual
        %function
        Ae_star(z,:)= getAstar( lambdastar , Ae_min, Ae_max,l_0_global(z,:),q0e(z,:), Lj(z,:),alpha_j(z,:) );
        
        %Compute element forces
        stress(z,:) = bar2s(ex(z,:),ey(z,:),[E Ae_star(z,:)],ed(z,:))/Ae_star(z,:);
        
        
    end
    
    %Compute the residual between current and previous iteration
    %regarding the updated area.
    res=Ae_star-A_global;
    Nres = norm(res);
    A_lastlast = A_last;
    A_last = A_global;
    A_global = Ae_star;
    
    
    %Plot vectors
    Ae_star_plot(:,iter) = Ae_star;
    Iter_plot(:,iter) = iter;
    g0(:,iter) = f'*a;
    stress_global(:,iter) = stress;
    res_plot(:,iter) = Nres;
    g1(:,iter) = A_global'*l_0_global-V_max;
    
    
    V = A_global'*l_0_global;
    
    fprintf('Iter %i Res %e\n',iter,Nres)
    
    
end

figure
plot(Iter_plot,res_plot,'LineWidth',1.5)
xlabel('Iteration')
ylabel('Residual')



%%
fac = 4*10^6;      %Magnification of Cross Section area.
magnfac = 1*10^-10; %Magnification of Deformed structure.
figure
myeldisp2(ex,ey,ed,[1 2 0],magnfac,A_global,fac);
grid on
xlabel('x,[m]')
ylabel('y,[m]')
title(['Structural optimization with MMA approximation'],'FontSize', 14)


figure
plot(Iter_plot,g0,'LineWidth',2)
hold on

plot(Iter_plot,g1,'LineWidth',2)'
grid on
xlabel('Iterationer','FontSize', 14)
ylabel('g_i','FontSize', 14)
legend('g0- objective function','g1- constraint function')
title('g_i functions with MMA approximation','FontSize', 14)
set(gca,'FontSize',14)
xlim([0 200])


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