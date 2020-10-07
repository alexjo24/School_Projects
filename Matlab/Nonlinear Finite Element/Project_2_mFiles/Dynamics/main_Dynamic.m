clc
clear all
close all
tic

%Material parameters
[ bc,edof,ex,ey,ez,loaddof,ndof,th,nen,nelm,f_int,a,ed,stress,~,~,ee_global] = data_dynamic();


E = 10e9;
v = 0.3;
ep = [E v];
rho  = 1700;
t=0;

%Selection of what time step length to compute.
prompt = ['Different time step lengths to select from:\n'...
    'A: delta_t = 2e-4 s \nB: delta_t = 2.5e-4 s \n'...
    'C: delta_t = 1e-4 s \nD: delta_t = 5e-4 s\n'...
    'What time step length: '];
A = 1;
B = 2;
C = 3;
D = 4;
selection = input(prompt)


if selection == 1
    % Time step length 1
    
    % Total running time
    time = 0.02;
    
    % Load steps
    n_LoadStep = 100;
    n_LoadStep_unload = 5;
    
    % Initiate linear time increment divided by the number of load steps. 
    delta_t = time/n_LoadStep;



elseif selection == 2 
    % Time step length 2

    %Total running time
    time = 0.05;

    %Load steps
    n_LoadStep = 200;
    n_LoadStep_unload = 4;

    %Initiate linear time increment divided by the number of load steps. 
    delta_t = time/n_LoadStep;

elseif selection == 3
    % Time step length 3

    %Total running time
    time = 0.015;

    %Load steps
    n_LoadStep = 150;
    n_LoadStep_unload = 10;

    %Initiate linear time increment divided by the number of load steps. 
    delta_t = time/n_LoadStep;

elseif selection == 4
    % Time step length 4

    %Total running time
    time = 0.07;

    %Load steps
    n_LoadStep = 140;
    n_LoadStep_unload = 2;

    %Initiate linear time increment divided by the number of load steps. 
    delta_t = time/n_LoadStep;
    
end

% Tolerance
TOL = 10^-6;

% Initiate
K = sparse(ndof,ndof);            %Stiffness matrix
M = sparse(ndof,ndof);            %Mass matrix
f = zeros(ndof,1);          %External force vector
P = zeros(ndof,1);

% Acceleration of nodal displacements.
a_dotdot = zeros(ndof,1);
a_dot = zeros(ndof,1);

% Initate gamme and beta constants, depending on stable region.
gamma  = 1/2;
beta = 1/4;

% Initate c - parameters
c1 = (1/(beta*(delta_t^2)));
c2 = (1/(beta*delta_t));
c3 = ((1-2*beta)/(2*beta));
c4 = (gamma/(beta*delta_t));
c5 = ((gamma-beta)/(beta));
c6 = ((delta_t*(gamma-2*beta))/(2*beta));


% Assemble the Mass matrix 
for k=1:nelm

ec = [ex(k,:) ; ey(k,:)];
[ Me ] = plan3gm( ec,th,rho );

indx = edof(k,2:end);
M(indx,indx) = M(indx,indx)+Me;
end



% Inital loading
P(loaddof) = 60e3;
P(125) = f_int(125);
P(126) = f_int(126);
delta_P = P/n_LoadStep_unload;



%%  %Initate damping parameters
d1 =0; %increase how much damping you want...
C = d1*M;


%Plane strain condition.
D_hooke = hooke(2,E,v);
D_hooke(:,3) = [];
D_hooke(3,:) = [];


for n = 1:n_LoadStep
    iter = 0;
    %Unloading according to a ramp 
    if n < n_LoadStep_unload
       
          P = P - delta_P;
          
    else
        P=zeros(ndof,1);
    end
        
    f = P;
  
    % Initiate iteration quantities
    a_dotdotStar = c1*a + c2*a_dot + c3*a_dotdot;
    a_dotStar    = c4*a + c5*a_dot + c6*a_dotdot;
    a_dotdotPrim = a_dotdotStar + d1*a_dotStar;
    a_dotdot_n = a_dotdot;
 
    % Predictor
    a_dotdot = ((1-gamma)/gamma)*a_dotdot;
    a = a + delta_t*a_dot+((delta_t^2 )/2)*((1-2*beta)*a_dotdot_n+...
                (2*beta*a_dotdot));
    
 
    % Effective residual.
    Geff =(c1+d1*c4)*M*a +f_int-f-M*a_dotdotPrim;
    
    disp(['Loadstep ' ,num2str(n),'----------------------------'])
    NRES = norm(Geff);
   

    while NRES > TOL || iter == 0
        
        iter = iter + 1;
        
       % Reset global matrices:
        K = sparse(ndof,ndof);
        f_int = zeros(ndof,1);
       
        %Element loop regarding tangent stiffness matrix
        for i=1:nelm

         ec = [ex(i,:) ; ey(i,:)];        
          
         %Compute tangent stiffness matrix
         [ Ke ] = plan3ge( ec,th,D_hooke,ed(i,:),stress(i,:)');
            
         %Assemble element matrices to global matrices:
         indx = edof(i,2:end);
         K(indx,indx) = K(indx,indx)+Ke;
         
        end
        
        %Compute effective stiffness matrix
        Keff = K + ((gamma*delta_t)/(beta*(delta_t^2)))*C+...
                    (1/(beta*delta_t^2))*M;        

           %Compute incremental displacments
           da = solveq(Keff,-Geff,bc);
           
           %Update the displacements
           a = a + da;
           ed = extract(edof,a);
           
           %Update the acceleration of the displacements
           a_dotdot = c1*a - a_dotdotStar;
           
           %Update the velocity of the displacements
           a_dot =    c4*a - a_dotStar;
           
           %Element loop of strain, stress and internal forces
           for j=1:nelm
               
               ec = [ex(j,:) ; ey(j,:)];

             % Compute the strain tensor ee and stress.
             [ ee,~ ] = plan3gs( ec,ed (j,:));
       
             %Store the strain for each element in a global strain matrix
             ee_global(j,:) = ee; 
             
             %St. Venant Kirchhoff constitutive law
             stress(j,:) = D_hooke*ee;
 
             % Compute the internal force vector
             [ ef ] = plan3gf( ec,th,ed(j,:),stress(j,:)');
           
             % Assemble internal element forces
              indx = edof(j,2:end);
              f_int(indx) = f_int(indx)+ef;

           end

           
           % Check effective residual.
            Geff =(c1+d1*c4)*M*a +f_int-f-M*a_dotdotPrim;

            res = Geff;
            res(bc(:,1)) = 0;
            NRES = norm(res);
            disp([num2str(iter),' Residual ',num2str(NRES)])

    end
    clf      
%     eldraw2(ex,ey,[1 4 0]);
    eldisp2(ex,ey,ed,[1 2 0],1);
    grid on
  	title(['Load step:', num2str(n)],'FontSize', 20)
    xlabel('m','FontSize', 14)
    ylabel('m','FontSize', 14)
    set(gca,'FontSize',14)
    drawnow

    

    
    %Update the time
    t = t + delta_t;
   
    
    %Plot vector regarding different ndof
      t_plot_temp(:,n) = t;
      f_plot2(:,n) = f(loaddof);
      f_plot125(:,n) = f(125);
      f_plot126(:,n) = f(126);


     %Compute kinetic and internal energy
     [KinE(:,n), IntE(:,n)] = plan3gEn( a_dot, M , D_hooke , ee_global',nelm, ex,ey, th );
      
     
     %Compute total energy
     TotalE(:,n) = KinE(:,n) + IntE(:,n);

end

%%

%Add reference step to plot
f_plot2 = [60000 ,f_plot2];
f_plot125 = [-60000, f_plot125];
f_plot126 = [-92047.0960688909  ,f_plot126];
t_plot = [0  ,t_plot_temp];


%Plot deformed and undeformed configuration of frame structure
figure
eldraw2(ex,ey,[1 4 0]);
eldisp2(ex,ey,ed,[1 2 0],1);
grid on


%Plot ramp curve in dof 2.
figure
plot(t_plot,f_plot2,'LineWidth',2)
grid on
xlabel('Time [s]','FontSize', 14)
ylabel('Force f [N]','FontSize', 14)
set(gca,'FontSize',14)
hold on

%Plot ramp curve in dof 125.
plot(t_plot,f_plot125,'LineWidth',2)
grid on
xlabel('Time [s]','FontSize', 14)
ylabel('Force f [N]','FontSize', 14)
set(gca,'FontSize',14)


%Plot ramp curve in dof 126.
plot(t_plot,f_plot126,'LineWidth',2)
grid on
xlabel('Time [s]','FontSize', 14)
ylabel('Force f [N]','FontSize', 14)
title('Force-time')
legend('1','125','126')
set(gca,'FontSize',14)
hold off


%Plot kinetic energy vs time
figure
subplot(2,1,1)
plot(t_plot_temp, KinE,'LineWidth',2)
grid on
title('Kinetic energy - Time','FontSize', 14)
set(gca,'FontSize',14)
hold on

%Plot internal energy vs time
plot(t_plot_temp,IntE,'LineWidth',2)
grid on
xlabel('Time [s]','FontSize', 14)
ylabel('Energy [J]','FontSize', 14)
title('Energy - Time','FontSize', 14)
set(gca,'FontSize',14)
legend('Kinetic energy','Internal energy')

%Plot total energy energy vs time
subplot(2,1,2)
plot(t_plot_temp,TotalE,'LineWidth',2)
grid on
xlabel('Time [s]','FontSize', 14)
ylabel('Total energy [J]','FontSize', 14)
title('Total energy  - Time','FontSize', 14)
set(gca,'FontSize',14)


toc


