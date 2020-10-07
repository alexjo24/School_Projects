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

%Material parameters
E = 200e09;                    %Young's modulus
rho_max = 1;              %Maximum density factor
rho_min = 10^-5;          %Minimum density factor
thickness = 20*10^-3;
V_box = 300*100*(10^-6)*thickness;
V_max = 0.4*V_box;        %Maximum allowed volume
v = 0.3;                  %Poissons tal


%Filter radius, according to Sigmund
%R = 0.0025*2;   %Fine mesh
R = 0.0071*2;   %Coarse mesh

%Initiate
K = sparse(ndof,ndof);
x_e = ones(nelm,1)*10^-1; %Introduce element Area vector, with all Area being
a_e = zeros(nelm,1);
p = 3;                %Penalization factor

ptype = 1; % Plane stress
ep =[ptype thickness 2];

%Hooke matrix
% D0 = (E/(1-v^2))*[1, v,    0;
%     v, 1,    0;
%     0, 0, (1-v)/2];

D0 = hooke(1,E,v);


%Initate Global cell for normalized stiffness matrix
Cell_K0e = cell(nelm,1);


for n=1:nelm
    
    a_e(n,:) = ElementArea( ec(n,:) ); %Area per element
    Vol_el(n,:) = thickness*ElementArea(ec(n,:));  %Volume per element
    
    %Initate constant norm element stiffness matrix
    [Cell_K0e{n}]=plani4e(exC(n,:),eyC(n,:),ep,D0);
    
    
    y_centrum=min(eyC(n,:))+(max(eyC(n,:))-min(eyC(n,:)))/2;
    x_centrum=min(exC(n,:))+(max(exC(n,:))-min(exC(n,:)))/2;
    kord_centrum=[x_centrum,y_centrum]';
    count=2;
    N(:,1)=edof(:,1);
    weights(:,1)=ones(max(size(edof)),1)*R;
    
    for lite=1:nelm
        
        y_centrum_test=min(eyC(lite,:))+(max(eyC(lite,:))-min(eyC(lite,:)))/2;
        x_centrum_test=min(exC(lite,:))+(max(exC(lite,:))-min(exC(lite,:)))/2;
        kord_test=[x_centrum_test,y_centrum_test]';
        
        test_val=norm(kord_test-kord_centrum);
        
        
        if test_val< R && lite ~= n
            weights(n,count)=R - test_val; % Cone-shaped function
            N(n,count)=lite;
            count=count+1;
        end
    end
    
    
    
end



for o = 1:nelm
    bot = 0;
    top = weights(o,1);
    
    for e = 2:min(size(weights))
        bot = bot + weights(o,e);
    
    end
    
    drhoTilde_drho(o,:) = top/bot; 
        
    
end


 

%Initate Tolerance
Tol = 10^-2;

%Initate quantaties
iter = 0;
Nres = 0;
alpha = 2; 


lambda_min=10^-20;
lambda_max=10^20;


%Optimization loop
while Nres > Tol || iter == 0
    
    iter = iter + 1;
    disp(['Iteration ' ,num2str(iter),'----------------------------'])
    
    K = sparse(ndof,ndof);
    
    %Element loop to assemble global stiffness matrix
    for i=1:nelm
        
        sum_top = 0;
        sum_bot = 0;
        
        for var=1:min(size(N))
            
            if N(i,var)~=0
                sum_top = sum_top+weights(i,var)*Vol_el(N(i,var))*x_e(N(i,var));
                
                sum_bot = sum_bot+weights(i,var)*Vol_el(N(i,var));
            end
        end
        
        x_e_tilde(i,:)=sum_top/sum_bot;
        
        D = (x_e(i,:)^p)*D0;
        
        Ke = plani4e(exC(i,:),eyC(i,:),ep,D);
        
        %Assemble element matrices to global matrices:
        indx = edof(i,2:end);
        K(indx,indx) = K(indx,indx)+Ke;
    end
    
    %Obtain displacement
    a = solveq(K,F,bc);
    ed = extract(edof,a);
    
    %     counter = 0;
    
    
    for j=1:nelm
        

        
        
        %Compute be according to eq 9.4
        
        var1 = (1/alpha)*(ed(j,:)*p*x_e_tilde(j,:)^(p-1)*Cell_K0e{j}*ed(j,:)'*(x_e(j,:)^(1+alpha)));
        
        be(j,:) = (1/alpha)*(ed(j,:)*p*x_e_tilde(j,:)^(p-1)*Cell_K0e{j}*drhoTilde_drho(j,:)*ed(j,:)'*(x_e(j,:)^(1+alpha)));
        
        
    end
    
    %Check derivative for each element fulfuill requirement, less than zero
    %         if counter == nelm
    %           prompt = 'All element dc_dxi are negative'
    %         else
    %           prompt = 'Some elements have dc_dxi positive'
    %           break
    %         end
    
    
    %         options = optimset('TolX',1e-12);
    %         lambdastar = fminbnd(@(lambda) phi_lambda...
    %         ( lambda,rho_max,rho_min,alpha,be,a_e,V_max,nelm )...
    %         ,lambda_min, lambda_max,options);
    
    
    
    lambdastar = fzero(@(lambda) dphidlambda_b...
        (lambda,rho_max,rho_min,alpha,be,a_e,V_max,nelm,thickness )...
        ,[lambda_min lambda_max]);
    
    for z=1:nelm
        
        %Find KKT point
        xe_star(z,:)= getx_e_star( lambdastar ,rho_max,rho_min,alpha,be(z,:),a_e(z,:) );
        
        
    end
    
    %      K=sparse(ndof,ndof);
    %Element loop to assemble global stiffness matrix
    %       for i=1:nelm
    %
    %           D = (((xe_star(i,:)^p)*E)/(1-v^2))*[1, v,    0;
    %                                           v, 1,    0;
    %                                           0, 0, (1-v)/2];
    %
    %           Ke = plani4e(exC(n,:),eyC(n,:),ep,D);
    %
    %          %Assemble element matrices to global matrices:
    %              indx = edof(i,2:end);
    %              K(indx,indx) = K(indx,indx)+Ke;
    %       end
    %
    %
    %
    %       %Obtain displacement
    %           a = solveq(K,F,bc);
    %           ed = extract(edof,a);
    
    %  for z=1:nelm
    %
    %         D = (((xe_star(z,:)^p)*E)/(1-v^2))*[1, v,    0;
    %                                             v, 1,    0;
    %                                             0, 0, (1-v)/2];
    %        %Compute element forces
    %       stress = plani4s(exC(z,:),eyC(z,:),ep,D,ed(z,:));
    %
    %       stress_global(z,:) = [stress(1,1) ,stress(2,2)  , stress(1,2) ];
    %
    %  end
    
    
    
    Nres = norm(xe_star-x_e);
    x_e = xe_star;
    
    %Plot vectors
    Iter_plot(:,iter) = iter;
    g0(:,iter) = F'*a;
    g1(:,iter) = x_e'*a_e*thickness-V_max;
    
    
    
    fprintf('Res %e\n',Nres)
    
    
    if iter == 200
        
        fprintf('Force break optimization loop at iteration: %i \n',iter)
        break
    end
end

%%

figure
fac = 5;       %Magnification of Cross Section area.
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
fill(exC',eyC',x_e_tilde)
colormap('jet')
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