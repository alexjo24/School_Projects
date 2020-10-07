clc
clear all
close all
tic

%Material parameters
[ bc,df_value,edof,ex,ey,ez,loaddof,ndof,th] = data();


E = 10e9;
v = 0.3;

nen = 3;              % Number of element nodes
nelm = max(edof(:,1)); 

%Load steps
n_LoadStep = 60;

%Tolerance
TOL = 10^-6;

%Initiate
K = sparse(ndof);            %Stiffness matrix
f = zeros(ndof,1);          %External force vector
f_int = zeros(ndof,1);      %Internal force vector
df = zeros(ndof,1);
stress = zeros(nelm,3);
ed = zeros(nelm,6);
u = zeros(ndof,1);

%Dof load increment, dof 2 is used.
df(loaddof)=(df_value/10);


%Plane strain condition.
D = hooke(2,E,v);
D(:,3) = []; 
D(3,:) = []; 


for n = 1:n_LoadStep
    
    iter = 0;
    
    %Add incremenet of external forces
    f = f + df;
    
    %Compute residual
    r = f - f_int;

    
    disp(['Loadstep ' ,num2str(n),'-----------------------------------'])
    NRES = norm(r);
     
    while NRES > TOL || iter == 0
        
        iter = iter + 1;
        
       %Reset global matrices:
       K = sparse(ndof,ndof);
       f_int = zeros(ndof,1);
       
 %% Bar element loop
        for i=1:nelm
 

          ec = [ex(i,:) ; ey(i,:)]';
           
           
          %Compute tangent stiffness matrix
          [ Ke ] = plan3ge( ec,th,D,ed(i,:),stress(i,:)');
           
         %Assemble element matrices to global matrices:
         indx = edof(i,2:end);
         K(indx,indx) = K(indx,indx)+Ke;
         
        end
        
        
           %Compute incremental displacments
           delta_u = solveq(K,r,bc);
           
           %Update the displacements
           u = u + delta_u;
           ed = extract(edof,u);
           
  %% Three node element loop, Green's strain and internal forces  
           for j=1:nelm
               
              ec = [ex(j,:) ; ey(j,:)];
                
             %Compute the strain tensor ee and stress.
             [ ee,~ ] = plan3gs( ec,ed(j,:));
                
             %Store strain in a global strain matrix.
             ee_glob(j,:)=ee';
             
             %Linear St. Venant-Kirchhoff const model , (2-D)
             stress(j,:) = D*ee;
             

             %Compute the internal force vector
             [ ef ] = plan3gf( ec,th,ed(j,:),stress(j,:)' );
           
             %Assemble internal element forces
              indx = edof(j,2:end);
              f_int(indx) = f_int(indx)+ef;
              
             
           end
           
          %Compute residual
           r = f - f_int;

        
            res= r;
            res(bc(:,1)) = 0;
            NRES = norm(res);
            fprintf('Iter %i Res %e\n',iter,NRES)
            

    end
   
    
end


figure
eldraw2(ex,ey,[1 4 0]);
xlabel('m','FontSize', 14)
ylabel('m','FontSize', 14)
set(gca,'FontSize',14)
title('Static loading')
grid on


figure
eldisp2(ex,ey,ed,[1 2 0],1);
xlabel('m','FontSize', 14)
ylabel('m','FontSize', 14)
set(gca,'FontSize',14)
title('Static loading')
grid on


toc


