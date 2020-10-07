clc
clear all
close all

[bc0,df_value,edof,edofB,ex,exB,ey,eyB,ez,ezB,loaddof,ndof,rc,th,xcen,ycen] = data();

%Material parameters
E = 10e9;
A_bar = 1;
v = 0.3;
ep = [E v];
epB = [E A_bar];
k = 1e7;


nelm =  max(edof(:,1));
nelmB = max(edofB(:,1));


%Load steps
n_LoadStep =50;

%Tolerance
epsilon = 10^-5;


%Initiate
K = sparse(ndof,ndof);            %Stiffness matrix
f = zeros(ndof,1);          %External force vector
f_int = zeros(ndof,1);      %Internal force vector
ed_tri = zeros(nelm,6);
ed_bar = zeros(nelmB,6);
u = zeros(ndof,1);
Stress = zeros(nelm,3);
eeB_global = zeros(nelmB,1);
forceB_global = zeros(nelmB,1);
df=zeros(ndof,1);
D_bar =zeros(nelmB,1);

%Plane strain condition.
D_hooke = hooke(2,E,v);
D_hooke(:,3) = [];
D_hooke(3,:) = [];
   
df(loaddof)=df_value;
 
    
for n = 1:n_LoadStep
    
    iter = 0;
    
    %Add incremenet of external forces
    f = f + df;

    %Compute residual
    r = f - f_int;
    
    bc = bc0;
    
    disp(['Loadstep ' ,num2str(n),'-----------'])
    NRES = norm(r);
    TOL = epsilon;
    
    while NRES > TOL || iter == 0
        
        iter = iter + 1;
        
       %Reset global matrices:
       K = sparse(ndof,ndof);
       f_int = zeros(ndof,1);
        
       
%% Three node element loop 
       for v=1:nelm

        ec = [ex(v,:) ; ey(v,:)]';
            
        %Compute element stiffness matrix for a three node element
       [Ke_tri] = plan3ge( ec,th,D_hooke,ed_tri(v,:),Stress(v,:)');
         
         % Assemble element matrices to global matrices
         indx = edof(v,2:end);
         K(indx,indx) = K(indx,indx)+Ke_tri;
         
       end
        
        
%% Bar element loop
        for i=1:nelmB
            ecB = [exB(i,:);
                   eyB(i,:);
                   ezB(i,:)];
         
               
          %Compute green's strain
         [eeB_global(i,:)] = bar3gs( ecB,ed_bar(i,:));
         
         %Calulate tangent stiffness tensor for bar element. 
         %Penalty method is implemented
         D_bar(i,:) = stiffb(ecB,eeB_global(i,:),k,rc);        
        
         ep = [D_bar(i,:) A_bar];  
        
         if D_bar(i,:)~=0
      
         %Compute element stiffness matrix for a bar element
         [Ke_bar] = bar3ge( ecB,ep,ed_bar(i,:),forceB_global(i,:)); 
         
         %Assemble element matrices to global matrices:
         indx = edofB(i,2:end);
         K(indx,indx) = K(indx,indx)+Ke_bar;
         else
             Ke_bar=zeros(6,6);
         end
         
         
         
        end
 
         
           %Compute incremental displacments
           delta_u = solveq(K,r,bc);
           
           %Update the displacements
           u = u + delta_u;
           
           
           
           
           ed_tri = extract(edof,u);
           ed_bar = extract(edofB,u);
           
          
           
           
  %% Three node element loop, Green's strain and internal forces          
           for j=1:nelm
               
              ec = [ex(j,:) ; ey(j,:)]';
  
             %Compute the strain tensor ee and stress.
             [ee_Global(j,:),eff ] = plan3gs( ec,ed_tri(j,:));
             
             %Linear St. Venant-Kirchhoff const model , (2-D)
             Stress(j,:) = D_hooke*ee_Global(j,:)';

             %Compute the internal force vector
             [ ef ] = plan3gf( ec,th,ed_tri(j,:),Stress(j,:) );
          
             %Assemble internal element forces
              indx = edof(j,2:end);
              f_int(indx) = f_int(indx)+ef;
 
           end
 
           
  %% Bar node element loop, Green's strain and internal forces
           for z=1:nelmB
           
               ecB = [exB(z,:);
                      eyB(z,:);
                      ezB(z,:)]; 

             %Compute the strain ee and normal force N.
             [eeB_global(z,:)] = bar3gs( ecB,ed_bar(z,:));
             
         
             %Compute normal force for each bar.
             %Penalty method is implemented
              [forceB_global(z,:)] = norfb(ecB,eeB_global(z,:),k,rc);    

            

             %Compute the internal force vector
             [ efB ] = bar3gf( ecB,ed_bar(z,:),forceB_global(z,:));
           
             %Assemble internal element forces
              indx = edofB(z,2:end);
              f_int(indx) = f_int(indx)+efB';
           end
           
            
           %Compute residual
            r = f - f_int;

            r(bc(:,1)) = 0;
            NRES = norm(r);
            disp([num2str(iter),' Residual ',num2str(NRES)])

            
            
    end
    
    
    
Force_bar=forceB_global';

clf
%Triangular element 
% eldraw2(ex,ey,[1 4 0]);
eldisp2(ex,ey,ed_tri,[1 2 0],1);
grid on
hold on

%Bar element
% eldraw2(exB,eyB,[1 4 0]);
eldisp2(exB,eyB,ed_bar,[1 5 0],1);
grid on


% Cylinder
Centers = [xcen, ycen];
viscircles(Centers, rc)
title(['Load step:', num2str(n)],'FontSize', 20)
xlabel('m','FontSize', 14)
ylabel('m','FontSize', 14)
set(gca,'FontSize',14)
drawnow


        
f_plot(:,n) = f(loaddof);
u_plot(:,n) = u(loaddof);


end

%Triangular element 
figure
eldraw2(ex,ey,[1 4 0]);
eldisp2(ex,ey,ed_tri,[1 2 0],1);
grid on
hold on

%Bar element
eldraw2(exB,eyB,[1 4 0]);
eldisp2(exB,eyB,ed_bar,[1 2 0],1);
grid on


% Cylinder
Centers = [xcen, ycen];
viscircles(Centers, rc)



%Plot internal force vs displacement
plot(u_plot,f_plot,'LineWidth',2)
grid on
xlabel('displacment [m]','FontSize', 14)
ylabel('Force [N]','FontSize', 14)
title('Force - Displacement','FontSize', 14)
set(gca,'FontSize',14)










