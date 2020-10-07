function [ bc,edof,ex,ey,ez,loaddof,ndof,th,nen,nelm,f_int,u,ed,stress,delta_u,f,ee_glob] = data_dynamic()
%Load topology

        load dyn_en.mat
        load ee_glob.mat
        bc=[2 0];

end

