clc
clear all
close all
format long
tic

load MBBFineMesh

nen = 4;

[ex,ey]=coordxtr(edof,coord,dof,nen);
ec = [ex,ey];

%Number of elements
nelm = max(edof(:,1));
ndof = max(max(dof));

%Material parameters
E = 200e9;                      % Young's modulus
rho_max = 1;                    % Maximum density factor
rho_min = 10^-5;                % Minimum density factor
thickness = 20*10^-3;           % Thickness
V_box = 300*100*(10^-6)*thickness;
V_max = 0.4*V_box;              % Maximum allowed volume
v = 0.3;                        % Poissons tal


%Filter radius, according to Sigmund
R = 0.0025*2.1;
%R = 0.0071*2;

%Initiate global matrix
K = sparse(ndof,ndof);          % Stiffness
K0=sparse(ndof,ndof);           % Normalized stiffness
Cell_K0e = cell(nelm,1);        % Normalized stiffness cell
x_e = ones(nelm,1)*10^-1;       % Initial design parameter
a_e = zeros(nelm,1);            % Area vector
p = 3;                          % Penalization factor

ptype = 1;                      % Plane stress
ep =[ptype thickness 2];

% Constitutive matrix
D0 = hooke(1,E,v);

for n=1:nelm
   
    % Calculate area and volume of each element
    a_e(n,:) = ElementArea( ec(n,:) );             
    Vol_el(n,:) = thickness*ElementArea(ec(n,:));  
    
    %Initate constant norm element stiffness matrix
    [Cell_K0e{n}]=plani4e(ex(n,:),ey(n,:),ep,D0);
    
   
    % Calculate centrum of current element and put in a vector
    y_centrum=min(ey(n,:))+(max(ey(n,:))-min(ey(n,:)))/2;
    x_centrum=min(ex(n,:))+(max(ex(n,:))-min(ex(n,:)))/2;
    kord_centrum=[x_centrum,y_centrum]';
    
    % Neigbour matrix 
    N(:,1)=edof(:,1);
    denom=0;
    
    for lite=1:nelm
        
        % Calculate centrum and check if neighbour 
        y_centrum_test=min(ey(lite,:))+(max(ey(lite,:))-min(ey(lite,:)))/2;
        x_centrum_test=min(ex(lite,:))+(max(ex(lite,:))-min(ex(lite,:)))/2;
        kord_test=[x_centrum_test,y_centrum_test]';
             
        if norm(kord_test-kord_centrum)< R %&& lite ~= n
            
            % Update neigbouring matrix 
            N(n,lite)=lite;     
            
            % Distance from centrum of current element cone-shaped function
            weights(n,lite)=R - norm(kord_test-kord_centrum);  
            
            % Global M-matrix 
            M(n,lite)=weights(n,lite)*Vol_el(n,:);
            denom=weights(n,lite)*Vol_el(n,:)+denom;
        end
    end
    M(n,:)=M(n,:)/denom;
end



%Initate Tolerance
Tol = 10^-2;

%Initate quantaties
iter = 0;
Nres = 0;
alpha = 2; %Conlin

% Interval to fzero 
lambda_min=10^-20;
lambda_max=10^20;

%Optimization loop
while Nres > Tol || iter == 0
    
    iter = iter + 1;
    disp(['Iteration ' ,num2str(iter),'----------------------------'])
    
    K = sparse(ndof,ndof);
    
    % Element loop to assemble global stiffness matrix
    for i=1:nelm
        
        x_e_tilde(i,:)=weights(i,:)*x_e/sum(weights(i,:));
        
%         sum_top = 0;
%         sum_bot = 0;
        
%         % Calculate the filtered design variable 
%         for hest=1:nelm
%             
%             if N(i,hest)~=0
%                 sum_top = sum_top+weights(i,hest)*Vol_el(N(i,hest))*x_e(N(i,hest));
%                 
%                 sum_bot = sum_bot+weights(i,hest)*Vol_el(N(i,hest));
%             end
%         end
%            

%         x_e_tilde(i,:)=sum_top/sum_bot;     
        
        % Constitutive matrix 
        D = (x_e_tilde(i,:)^p)*D0;
        
        % Stiffness matrix 
        Ke = plani4e(ex(i,:),ey(i,:),ep,D);
        
        %Assemble element matrices to global matrices:
        indx = edof(i,2:end);
        K(indx,indx) = K(indx,indx)+Ke;
    end
    
    %Obtain displacement
    a = solveq(K,F,bc);
    ed = extract(edof,a);
    
    % Calculate the sensitivity 
    
    
    % Assemble and calculate constants to the approximation 
    for j=1:nelm
        dC_drhotilde(j,:)= ed(j,:)*(-p*x_e_tilde(j,:)^(p-1))*Cell_K0e{j}*ed(j,:)';                
    end
    
    dC_drho=dC_drhotilde'*M;
    
 
    be= -(1/alpha)*x_e.^(1+alpha).*dC_drho';
  
    % Solve the dual Lagrangian function 
    lambdastar = fzero(@(lambda) dphidlambda_b...
        (lambda,rho_max,rho_min,alpha,be,a_e,V_max,nelm,thickness )...
        ,[lambda_min lambda_max]);
    
    % Find KKT point
    for z=1:nelm     
        xe_star(z,:)= getx_e_star( lambdastar ,rho_max,rho_min,alpha,be(z,:),a_e(z,:) );
    end  
    
    Nres = norm(xe_star-x_e);
    x_e = xe_star;
    
    %Plot vectors
    Iter_plot(:,iter) = iter;
    g0(:,iter) = F'*a;
    g1(:,iter) = x_e_tilde'*a_e*thickness-V_max;
    
    fprintf('Res %e\n',Nres)

    if iter == 30
        fprintf('Force break optimization loop at iteration: %i \n',iter)
        break
    end
end

%%

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
colormap('jet')
fill(ex',ey',x_e_tilde)
title('Density Filter')
xlabel('x [m]')
ylabel('y [m]')

toc