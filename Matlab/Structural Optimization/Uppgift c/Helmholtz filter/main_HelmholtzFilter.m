clc
clear all
close all
clear vars
format long


el_length_Fine = 0.0025;
el_length_Coarse = 0.0071;

%Selection of coarse mesh or fine mesh.
prompt = ['Select Coarse Mesh or Fine Mesh:\n'...
    'A1: Coarse Mesh \nB1: Fine Mesh \n'];
A1 = 1;
B1 = 2;
selection = input(prompt)


if selection == 1
    
    load MBBCoarseMesh
    
    %Selection of what length scale R to use.
    prompt = ['Select a length scale R\n'...
        'A2: R = Element Side*1  \nB2: R = Element Side*2 \n'...
        'C2: R = Element Side*2.5\n'];
    A2 = 1;
    B2 = 2;
    C2 = 3;
    selection2 = input(prompt)
    
    
    
    if selection2 == 1
        
        R = el_length_Coarse;
        
    elseif selection2 == 2
        
        R = el_length_Coarse*2;
        
    elseif selection2 == 3
        
        R = el_length_Coarse*2.5;
        
        
    end
    
elseif selection == 2
    
    load MBBFineMesh
    
    %Selection of what length scale R to use.
    prompt = ['Select a length scale R\n'...
        'A2: R = Element Side*1  \nB2: R = Element Side*2 \n'...
        'C2: R = Element Side*2.5\n'];
    A2 = 1;
    B2 = 2;
    C2 = 3;
    selection2 = input(prompt)
    
    
    
    if selection2 == 1
        
        R = el_length_Fine;
        
    elseif selection2 == 2
        
        R = el_length_Fine*2;
        
    elseif selection2 == 3
        
        R = el_length_Fine*2.5;
        
        
    end
    
end


%Input min and max values for fzero function.
lambda_min=10^-8;
lambda_max=10^-1;

nen = 4; %Number of nodes per element

%Extract x and y coordinates
[exC,eyC]=coordxtr(edof,coord,dof,nen);
ec = [exC,eyC];


nelm = max(edof(:,1));  %Number of elements
ndof = max(max(dof));   %Number of degree of freedom, displacement part
nnod = ndof/2;          %Number of nodes

%Material parameters
E = 200e9;               %Young's modulus
rho_max = 1;              %Maximum density factor
rho_min = 10^-5;          %Minimum density factor
thickness = 20*10^-3;     %Thickness
V_box = 300*100*(10^-6)*thickness;
V_max = 0.4*V_box;        %Maximum allowed volume
v = 0.3;                  %Poissons tal
eq = 1;


%Initiate
x_e = ones(nelm,1)*10^-3; %Introduce element Area vector, with all Area being
a_e = zeros(nelm,1);      %Introduce Area vector
T = sparse(nnod,nelm);    %Introduce T matrix, [nnod x nelm]
M = sparse(nnod,nnod);    %Introduce Mass matrix
K_dens = sparse(nnod,nnod); %Introduce constant stiffness matrix,Ke [4x4]
K_rho =sparse(ndof,ndof); %Introduce stiffness matrix, Ke [8x8]
dgdrhotilde = sparse(nnod,1);
p = 3;                    %Penalization factor

ptype = 1; %Plane stress
ep =[ptype thickness 2];

%Hooke matrix
D0 = hooke(ptype,E,v);


%Parameters for flw2i4e
ep2 = [thickness 2];
D2 = [1 0;
    0 1];


for n=1:nelm
    
    a_e(n,:) = ElementArea( ec(n,:) );   %Area per element
    Vol_el(n,:) = thickness*ElementArea(ec(n,:));  
    
    
    %Compute T matrix
    [Ke,Te] = flw2i4e(exC(n,:),eyC(n,:),ep2,D2,eq);
    
    %Assemble element matrices to global matrices:
    indx = enod(n,2:end);
    T(indx,n) = T(indx,n)+Te;
    
    indx = enod(n,2:end);
    K_dens(indx,indx) = K_dens(indx,indx)+Ke.*R^2;
    
    
    %Compute Mass matrix
    Me = flw2i4m(exC(n,:),eyC(n,:),thickness);
    
    %Assemble element matrices to global matrices:
    indx = enod(n,2:end);
    M(indx,indx) = M(indx,indx)+Me;
    
    
    
end

%Compute constant drhoTilde_drho matrix.
drhoTilde_drho = (K_dens + M)\T;


%Initate Tolerance
Tol = 10^-9;

%Initate quantaties
iter = 0;
Nres = 0;
alpha = 2;


%Optimization loop
while Nres > Tol || iter == 0
    
    iter = iter + 1;
    disp(['Iteration ' ,num2str(iter),'----------------------------'])
    
    %Reset global matrices
    K_rho = sparse(ndof,ndof);
    dgdrhotilde = sparse(nnod,1);
    
    %Compute filtered density
    x_e_tilde = drhoTilde_drho*x_e;
    
    %Extract the densities at the nodal points.
    ed_rho=extract(enod,x_e_tilde);
    
    
    
    %Element loop to assemble global stiffness matrix
    for i=1:nelm
        
        Ke_rho=plani4e_rho(exC(i,:),eyC(i,:),ep,D0,ed_rho(i,:),p);
        
        %Assemble element matrices to global matrices:
        indx = edof(i,2:end);
        K_rho(indx,indx) = K_rho(indx,indx)+Ke_rho;
        
        
    end
    
    %Obtain displacement
    a = solveq(K_rho,F,bc);
    ed = extract(edof,a);
    
    
    for j=1:nelm
        
        %Compute dgdrhotilde
        dgdrhotilde_e=getdgdrhotilde_el(exC(j,:),eyC(j,:),ep,D0,ed(j,:)...
            ,ed_rho(j,:),p);
        
        %Assemble element matrices to global matrices:
        indx = enod(j,2:end);
        dgdrhotilde(indx,1) = dgdrhotilde(indx,1)+dgdrhotilde_e;
        
        
    end
    
    %Compute global dg_drho
    dg_drho = dgdrhotilde'*drhoTilde_drho;
   
    
    b = (1/alpha).*dg_drho'.*(x_e.^(1+alpha));
    
    
    %Solve the dual function
    lambdastar = fzero(@(lambda) dphidlambda_b...
        (lambda,rho_max,rho_min,alpha,b,a_e,V_max,nelm,thickness )...
        ,[lambda_min lambda_max]);
    
    for z=1:nelm
        
        %Find new rho value per element based on new lambda from the dual
        %function
        xe_star(z,:)= getx_e_star( lambdastar ,rho_max,rho_min,alpha,b(z,:),a_e(z,:) );
        
    end
    
    
    Nres = norm(xe_star-x_e);
    x_e = xe_star;
    
    %Plot vectors
    Iter_plot(:,iter) = iter;
    g0(:,iter) = F'*a;
    g1(:,iter) = x_e'*a_e*thickness-V_max;
    
    
    fprintf('Res %e\n',Nres)
    
    
    if iter == 1000
        
        fprintf('Force break optimization loop at iteration: %i \n',iter)
        break
    end
    
    
    if any(dg_drho<0) == 0
       fprintf('R to small')
    end
        
        
end

%%

ed_el = extract(enod,x_e_tilde);
figure
fill(exC',eyC',ed_el')
title('Helmholtz density PDE filter')
xlabel('x, [m]')
ylabel('y, [m]')
colormap('Jet')
axis([0 0.3 0 0.12])

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

