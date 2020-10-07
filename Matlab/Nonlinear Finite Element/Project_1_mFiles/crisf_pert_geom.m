clear all
close all
clc
tic

%Geometry
[ex,ey,ez,edof,ndof,nelm,plotdof,plotdof2,ep,P,bc0]=geom2016_pert();


 E=1;
 A_0=1;

v=0.3;
G_s=E/(2*(1+v));

n_load_steps=100;
TOL=10^-8;

%Initiate
K=zeros(ndof);
f = zeros(ndof,1); %External force vector
f_int = zeros(ndof,1); %Internal force vector
ed=zeros(nelm,6);
a_i = zeros(ndof,1);
a_0 = zeros(ndof,1);
lambda_i = 0;
lambda_0 = 0;
kraft=zeros(nelm,1);
stress = zeros(nelm,1);
eps=zeros(nelm,1);
f_plot_zero = 0;
u_plot_zero = 0;

%Load inc
%Introduce psi and l
l = 0.1;
psi = 1;


%Introduce perturbations
for n=1:n_load_steps
    iter=0;
    
    bc=bc0;
    f = lambda_0*P;
    

    G = f - f_int;
    nres=norm(G);
    disp([num2str(n),'------------------------'])
    
    
%     delta_a_i = zeros(ndof,1);
%     delta_lambda_i = 0;

    while nres>TOL || iter==0
        iter=iter+1;
        
        %Nollställa globala matriser
        K = zeros(ndof); %Stiffness matrix
        f_int = zeros(ndof,1);
      
        for j=1 : nelm         
            ec=[ex(j,:)
                ey(j,:)
                ez(j,:)]; 
            
           
            
            D= G_s*((1/(sqrt(2*eps(j,:)+1))) + (2/((2*eps(j,:)+1)^2)));
            ep_temp = [D, A_0];
            Ke = bar3ge_k(ec,ep_temp,ed(j,:),kraft(j,:));
            indx = edof(j,2:end);
            K(indx,indx)=K(indx,indx)+Ke;          
        end

        daG = solveq(K,G,bc);
        daP = solveq(K,P,bc);
     
        if iter == 1
                a1 = (daP'*daP) + psi*(P'*P);
                a2 = 0; 
                a3 = -l^2;
             
              if n == 1
  

                %Kan vara positiv eller negativ
                dlambda = sqrt(-a3/a1);
  
                else 
                
                delta_a_n = a_0 - a_0_1;   

                s=sign(delta_a_n'*daP);
                dlambda = s*l/(sqrt(daP'*daP+psi*P'*P));
            
                end
            else
             
            delta_a_i = a_i - a_0;
         
            delta_lambda_i = lambda_i-lambda_0;
         
            a1 = daP'*daP + psi*(P'*P);
            a2 = 2*daP'*(delta_a_i+daG)+2*psi*delta_lambda_i*(P'*P);
            a3 = (delta_a_i+ daG)'*(delta_a_i+daG)+psi*(delta_lambda_i^2)*(P'*P)-l^2;

            dlambda = roots([a1,a2,a3]);
            
            
            if isreal(dlambda) == 0
            
               %Change l-value manually
                error('Decrease load step length l, imaginary solution')

            end
            

            a4 = delta_a_i'*(delta_a_i+daG);
            a5 = delta_a_i'*daP;
            
             angle1 = acos((a4 + a5*dlambda(1))/l^2);
             angle2 = acos((a4 + a5*dlambda(2))/l^2);
             
             
                     if angle1 < angle2
                         dlambda = dlambda(1);
                     else
                         dlambda = dlambda(2);
                     end

            

        end
            
            
        
            lambda_i =  lambda_i + dlambda;
      
            a_i = a_i + daG + dlambda*daP;

            ed = extract(edof,a_i);
        
        for j=1: nelm
            ec=[ex(j,:)
                ey(j,:)
                ez(j,:)];
            
         
 D = G_s*((1/(sqrt(2*eps(j,:)+1))) + (2/((2*eps(j,:)+1)^2)));
            

              ep_temp = [D, A_0];
            [~,eps(j,:)] =bar3gs(ec,ep_temp,ed(j,:));

             stress(j,:)=stress1D(E,v,ec,ed(j,:));   
      
              kraft(j,:)=A_0*stress(j,:);     

            ef = bar3gf( ec, ed(j,:),kraft(j,:));
            
            indx = edof(j,2:end);
            f_int(indx)=f_int(indx)+ef';
       
        end 
    
        f = lambda_i*P;
        G = f - f_int;
      
        res = G;
        
        res(bc(:,1)) = 0;
        nres=norm(res);
        disp([num2str(iter),'  res ',num2str(nres)])
   
    end
       
    a_0_1 = a_0;
    a_0 = a_i;
           
    lambda_0 = lambda_i;
    
    
          f_plot1(n) = f(plotdof);
          f_plot2(n) = f(plotdof);
          u_plot1(n) = a_i(plotdof);
          u_plot2(n) = a_i(plotdof2);
          
          
          
           clf
           eldisp3(ex,ey,ez,ed,[1 4 1],1);
           title(['Load step:', num2str(n)],'FontSize', 20)
           grid on
           drawnow 
      
          
end


u_plot1 = [u_plot_zero, u_plot1];
u_plot2 = [u_plot_zero, u_plot2];
f_plot1 = [f_plot_zero, f_plot1];
f_plot2 = [f_plot_zero, f_plot2];

figure
plot(u_plot1,f_plot1);
title('Force-displacement, B-direction')
xlabel('Displacement [mm]')
ylabel('Force [N]')
grid on

figure 
plot(u_plot2,f_plot2);
title('Force-displacement, A-direction')
xlabel('Displacement [mm]')
ylabel('Force [N]')
grid on

figure 
plot(u_plot1,u_plot2);
title( 'A-direction and B-direction')
xlabel('Displacement [mm]')
ylabel('Displacement [mm]')
grid on

figure
eldraw3(ex,ey,ez,[1 2 1])
eldisp3(ex,ey,ez,ed,[1 4 1],1);
grid on


toc