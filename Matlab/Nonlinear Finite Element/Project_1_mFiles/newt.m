clear all
close all
clc
tic

%Geometry
[ex,ey,ez,edof,ndof,nelm,plotdof,plotdof2,ep,P,bc0]=geom2016();
E=1;
A_0=1;
v=0.3;

n_load_steps=40;
%delta_f=1*10^-4; 
TOL=10^-8;

K=zeros(ndof);
f = zeros(ndof,1); %External force vector
f_int = zeros(ndof,1); %Internal force vector

es=0;
kraft=zeros(nelm,1);
span=zeros(nelm,1);
eps=zeros(nelm,1);

ed=zeros(nelm,6);
u=zeros(ndof,1);


f_plot=zeros(n_load_steps,1);
u_plot_4=zeros(n_load_steps,1);
u_plot_5=zeros(n_load_steps,1);
u_plot_6=zeros(n_load_steps,1);
f_eff_plot = zeros(n_load_steps,1);
f_plot_zero = 0;
u_plot_zero = 0;
delta_f=P;

for n=1:n_load_steps
    iter=0;
%   f_n=f_(n-1)+delta_f_n
    bc=bc0;
    
    f=f+delta_f;
    
    res=f-f_int;
    nres=norm(res);
    disp([num2str(n),'------------------------'])
    
   % ed=extract(edof,u);
    while nres>TOL || iter==0
        iter=iter+1;
        
        %Nollställa globala matriser
        K = zeros(ndof); %Stiffness matrix
        f_int = zeros(ndof,1);
        
        for j=1 : nelm         
            ec=[ex(j,:)
                ey(j,:)
                ez(j,:)]; 
            Ke = bar3ge(ec,ep,ed(j,:),kraft(j,:));
            indx = edof(j,2:end);
            K(indx,indx)=K(indx,indx)+Ke;          
        end
        
        delta_u=solveq(K,res,bc);
    
        u=u+delta_u;

        ed=extract(edof,u);
        
        for j=1: nelm
            ec=[ex(j,:)
                ey(j,:)
                ez(j,:)];
        
            [~,eps(j,:)] =bar3gs(ec,ep,ed(j,:));
            span(j,:)=eps(j,:)*E;           
            kraft(j,:)=A_0*span(j,:);

               
            ef = bar3gf( ec, ed(j,:),kraft(j,:));
            
            indx = edof(j,2:end);
            f_int(indx)=f_int(indx)+ef';
       
        end 
    
        res=f-f_int;
        res(bc(:,1)) = 0;
        nres=norm(res);
        disp([num2str(iter),'  res ',num2str(nres)])
   
    end
       
          f_plot1(n) = f(plotdof);
          f_plot2(n) = f(plotdof2);
          u_plot1(n) = u(plotdof);
          u_plot2(n) = u(plotdof2);

end

u_plot1 = [u_plot_zero, u_plot1];
u_plot2 = [u_plot_zero, u_plot2];
f_plot1 = [f_plot_zero, f_plot1];
f_plot2 = [f_plot_zero, f_plot2];



figure
plot(u_plot1,f_plot1,'LineWidth',2);

hold on

 
plot(u_plot2,f_plot2,'LineWidth',2);
title('Force-displacement')
xlabel('Displacement [mm]')
ylabel('Force [N]')
legend('B-direction', 'A-direction')
grid on

figure
eldraw3(ex,ey,ez,[1 2 1])
eldisp3(ex,ey,ez,ed,[1 4 1],1);
grid on


toc