function [Ex,Ey,Ez,edof,ndof,nelm,plotdof,plotdof2,ep,P,bc]=geom2016_pert()

% Coordinates
x1=1.697;
x2=1.414;

coord=[-x1 -1 0  %1 node number
        0  -1 0  %2
        x1 -1 0  %3
       -x1  1 0  %4
        0   1 0  %5
        x1  1 0  %6
       -x1  0 1  %7
        0   0 1  %8
        x1  0 1+0.5];%9   

% Degrees of freedom for each node
dof=[1 2 3      %1 node number
     4 5 6      %2
     7 8 9      %3   
     10 11 12   %4
     13 14 15   %5
     16 17 18   %6
     19 20 21   %7
     22 23 24   %8
     25 26 27]; %9

% Topology matrix
edof=[1  dof(1,:) dof(7,:)
      2  dof(2,:) dof(7,:)
      3  dof(2,:) dof(8,:)
      4  dof(2,:) dof(9,:)
      5  dof(3,:) dof(9,:)
      6  dof(4,:) dof(7,:)
      7  dof(5,:) dof(7,:)
      8  dof(5,:) dof(8,:)
      9  dof(5,:) dof(9,:)
      10 dof(6,:) dof(9,:)
      11 dof(7,:) dof(8,:)
      12 dof(8,:) dof(9,:)];

% Extract element coordinates
[Ex,Ey,Ez]=coordxtr(edof,coord,dof,2);

% figure(1);
eldraw3(Ex,Ey,Ez,[1 5 1]);


% Essential boundary conditions
bc=[dof(1,:)' [0 0 0]'
    dof(2,:)' [0 0 0]'
    dof(3,:)' [0 0 0]'
    dof(4,:)' [0 0 0]'
    dof(5,:)' [0 0 0]'
    dof(6,:)' [0 0 0]'];

% pause
% global vectors sizes
ndof=max(max(dof));
nelm=max(edof(:,1));

% External load increment
P = zeros(ndof,1); 
P([dof(7,3) dof(8,3) dof(9,3)])=-0.03*[1.5 1 1.5];

plotdof=dof(8,3);
plotdof2=dof(9,3);

Edof=edof;


% truss properties
ep=[1 1];
end