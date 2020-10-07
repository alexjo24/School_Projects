clc
clear all
close all
format long

load MBBCoarseMesh

nen = 4;

[exC,eyC]=coordxtr(edof,coord,dof,nen);
ec = [exC,eyC];

%Number of elements
nelm = max(edof(:,1));
ndof = max(max(dof));
nnod = ndof/2;

%Material parameters
E = 200e09;                    %Young's modulus
rho_max = 1;              %Maximum density factor
rho_min = 10^-4;          %Minimum density factor
thickness = 20*10^-3;
V_box = 300*100*(10^-6)*thickness;
V_max = 0.4*V_box;        %Maximum allowed volume
v = 0.3;                  %Poissons tal
eq = 1;

%Filter radius, according to Sigmund
%R = 0.0025*2;   %Fine mesh
R = 0.0071*1.5;   %Coarse mesh

%Initiate
K = sparse(ndof,ndof);
x_e = ones(nelm,1)*10^-3; %Introduce element Area vector, with all Area being
a_e = zeros(nelm,1);
T = sparse(nnod , nelm);
M = sparse(nnod,nnod);
K_dens = sparse(nnod,nnod);
K_rho =sparse(ndof,ndof);
p = 3;                %Penalization factor

ptype = 1; % Plane stress
ep =[ptype thickness 2];

%Hooke matrix
D0 = (E/(1-v^2))*[1, v,    0;
                  v, 1,    0;
                  0, 0, (1-v)/2];


%Parameters for flw2i4e
ep2 = [thickness 2];
D2 = [1 0;
     0 1];

%Initate Global cell for normalized stiffness matrix
Cell_K0e = cell(nelm,1);


for n=1:nelm
    
    a_e(n,:) = ElementArea( ec(n,:) ); %Area per element
    Vol_el(n,:) = thickness*ElementArea(ec(n,:));  %Volume per element
    
    %Initate constant norm element stiffness matrix
  %  [Cell_K0e{n}]=plani4e(exC(n,:),eyC(n,:),ep,D0);
    
%     %Find neighbour elements to a center element and compute weights.
%     y_centrum=min(eyC(n,:))+(max(eyC(n,:))-min(eyC(n,:)))/2;
%     x_centrum=min(exC(n,:))+(max(exC(n,:))-min(exC(n,:)))/2;
%     kord_centrum=[x_centrum,y_centrum]';
%     count=2;
%     N(:,1)=edof(:,1);
%     weights(:,1)=ones(max(size(edof)),1)*R;
%     
%     for lite=1:nelm
%         
%         y_centrum_test=min(eyC(lite,:))+(max(eyC(lite,:))-min(eyC(lite,:)))/2;
%         x_centrum_test=min(exC(lite,:))+(max(exC(lite,:))-min(exC(lite,:)))/2;
%         kord_test=[x_centrum_test,y_centrum_test]';
%         
%         test_val=norm(kord_test-kord_centrum);
%         
%         
%         if test_val< R && lite ~= n
%             weights(n,count)=R - test_val; % Cone-shaped function
%             N(n,count)=lite;
%             count=count+1;
%         end
%     end
    
    %Compute T matrix
    [Ke,Te] = flw2i4e(exC(n,:),eyC(n,:),ep2,D2,eq);

    %Assemble element matrices to global matrices:
        indx = enod(n,2:end);
        T(indx,n) = T(indx,n)+Te;
    
        indx = enod(n,2:end);
        K_dens(indx,indx) = K_dens(indx,indx)+Ke*R^2;
        
        
        
    %Compute Mass matrix
     Me = flw2i4m(exC(n,:),eyC(n,:),thickness);
    
     indx = enod(n,2:end);
     M(indx,indx) = M(indx,indx)+Me;
    
     %drhoTilde_drho_cell{n}=(Ke + Me)\Te;

end
%drhoTilde_drho = ((K_dens + M)^-1)*T;
drhoTilde_drho=solveq(K_dens+M,T);

%Initate Tolerance
Tol = 10^-3;

%Initate quantaties
iter = 0;
Nres = 0;
alpha = 2; %Conlin


lambda_min=10^-20;
lambda_max=10^20;


%Optimization loop
while Nres > Tol || iter == 0
    
    iter = iter + 1;
    disp(['Iteration ' ,num2str(iter),'----------------------------'])
    
    %K = sparse(ndof,ndof);
    K_rho = sparse(ndof,ndof);
    
    x_e_tilde = drhoTilde_drho*x_e;
    ed_rho=extract(enod,x_e_tilde);
    
    %Element loop to assemble global stiffness matrix
    for i=1:nelm
        
%         sum_top = 0;
%         sum_bot = 0;
%         
%         for hest=1:min(size(N))
%             
%             if N(i,hest)~=0
%                 sum_top = sum_top+weights(i,hest)*Vol_el(N(i,hest))*x_e(N(i,hest));
%                 
%                 sum_bot = sum_bot+weights(i,hest)*Vol_el(N(i,hest));
%             end
%         end
        
%         x_e_tilde(i,:)=sum_top/sum_bot;
        
        D = (x_e(i,:)^p)*D0;
        
%         Ke = plani4e(exC(i,:),eyC(i,:),ep,D);
%         
%         %Assemble element matrices to global matrices:
%         indx = edof(i,2:end);
%         K(indx,indx) = K(indx,indx)+Ke;
        
        Ke_rho=plani4e_rho(exC(i,:),eyC(i,:),ep,D,ed_rho(i,:),p);
        
        %Assemble element matrices to global matrices:
        indx = edof(i,2:end);
        K_rho(indx,indx) = K_rho(indx,indx)+Ke_rho;
        
        
    end
    
    
    %Obtain displacement
    %a = solveq(K,F,bc);
    %ed = extract(edof,a);
    a=solveq(K_rho,F,bc);
    ed=extract(edof,a);
    
    
        dgdrhotilde=zeros(nnod,1);
   
        for z= 1:nelm
           %Compute dgdrhotilde
       
           dgdrhotilde_e=getdgdrhotilde_el(exC(z,:),eyC(z,:),ep,D,ed(z,:),ed_rho(z,:),p);
         %  dg_drho = dg_drho + dgdrhotilde'*drhoTilde_drho_cell{z}; 
            
            indx = enod(z,2:end);
            dgdrhotilde(indx) = dgdrhotilde(indx)+dgdrhotilde_e;

        end

        

       dc_drho=dgdrhotilde'*drhoTilde_drho;
       
        for j=1:nelm
         
            be(j,:) = (1/alpha)*(x_e(j,:)^(1+alpha))*dc_drho(:,j);
         
        end

    
    lambdastar = fzero(@(lambda) dphidlambda_b...
        (lambda,rho_max,rho_min,alpha,be,a_e,V_max,nelm,thickness )...
        ,[lambda_min lambda_max]);
    
    for z=1:nelm
        
        %Find KKT point
        xe_star(z,:)= getx_e_star( lambdastar ,rho_max,rho_min,alpha,be(z,:),a_e(z,:) );
        
        
    end
  
    
    Nres = norm(xe_star-x_e);
    x_e = xe_star;
    
    %Plot vectors
    Iter_plot(:,iter) = iter;
    g0(:,iter) = F'*a;
    g1(:,iter) = x_e'*a_e*thickness-V_max;
    
    
    
    fprintf('Res %e\n',Nres)
    
    
    if iter == 300
       
        
        fprintf('Force break optimization loop at iteration: %i \n',iter)
        break
    end
end

%%
figure
fac = 5;          %Magnification of Cross Section area.
magnfac = 10^-12; %Magnification of Deformed structure.
myeldisp2(exC,eyC,ed,[1 2 0],magnfac,x_e_tilde,fac);
grid on


figure
plot(Iter_plot,g0,'LineWidth',2)
grid on
xlabel('Iterationer')
ylabel('g_0')
title('Compliance history')


figure
plot(Iter_plot,g1,'LineWidth',2)
xlabel('Iterationer')
ylabel('g_1')
grid on
title('Weight history')

figure
ed_el = extract(enod,x_e_tilde);
fill(exC',eyC',ed_el')

%Check Area fulfills Box constraint

% A_max_check = max(x_e);
% A_min_check = min(x_e);
%
%     if A_max_check <= Ae_max
%         disp('Valid Box constraint regarding A_max')
%     else
%         disp('Box constraint not fulfilled')
%     end
%
%
%      if A_min_check >= Ae_min
%         disp('Valid Box constraint regarding A_min')
%         else
%         disp('Box constraint not fulfilled')
%      end
%
%
%   procent_min=length(find(x_e == Ae_min))/nelm
%   procent_max=length(find(x_e == Ae_max))/nelm