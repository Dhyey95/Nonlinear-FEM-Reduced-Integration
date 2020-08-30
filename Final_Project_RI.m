clc
clear all;

%Material Properties
E = 30000; %psi
v = 0.4999; %poisson's ratio
lamda = v*E/((1+v)*(1-2*v));
mu = E/(2*(1+v));
kron_delta = [1 0 0; 0 1 0; 0 0 1];
Csmat = Materialconst(lamda,mu,kron_delta);

%Mesh Generation
l_x = 50; %inch
l_y = 5; %inch
n_x = 81;
n_y = 9;
center_node = 0.5*n_x*(n_y+1);
ele_x = n_x - 1;
ele_y = n_y - 1;
nn = n_x*n_y;
nele = ele_x*ele_y;
x_coord = linspace(0,l_x,n_x);
y_coord = linspace(0, l_y, n_y);
n_coord = zeros(nn,2);
k = 1;
for j = 1:size(y_coord,2)
    for i = 1:size(x_coord,2)
        n_coord(k,1) = x_coord(i);
        n_coord(k,2) = y_coord(j);
        k = k+1;      	
    end
end

ele_number = zeros(nele, 4); %2D elements with 4 node in each element
i=1; 
for j=1:ele_y
    for k = j+((j-1)*(n_x -1)):(j+((j-1)*(n_x-1)))+(n_x-2)
        ele_number(i,:) = [k,k+1,k+n_x+1,k+n_x];
        i=i+1;
    end
end

%Gaussian Integration
% 1 point & 2 point gauss integration
Ip1 = 0; w1 = 2; Np1 = 1;
Ip2 = [-1/sqrt(3), 1/sqrt(3)]; w2 = [1 1]; Np2 = 2;

% Load applied
P = -250/l_y;
N_load_pt = 500;
tolerance = 10^-5;
delta_defl = zeros(2*nn,1);
dof_inds = @(x) [2*x-1; 2*x];
delta_P = P/N_load_pt;

for i = 1:N_load_pt 
    
    F_global_external = zeros(2*nn,1);
    
    Pnew =  i*delta_P;
     
    for j = 1:nele
        ele_dof = ele_number(j,:);
        ele_coord = n_coord(ele_number(j,:),:);
        x = ele_coord(:,1);
        y = ele_coord(:,2);
        
        Fgamma = zeros(8,1);
        if mod(j,(n_x-1))==0
            for k = 1:Np1
                eps = 1.0;
                eta = Ip1(k);
                [N,J,Bun] = shape(eps,eta,x,y);
                h_g = [0;Pnew];
                j_g = norm([J(2,2),J(2,1)]);
                Fgamma = Fgamma + (N')*h_g*j_g*2;
            end
        end
    
        
        global_inds = reshape(dof_inds(ele_dof),1,[]); 
        inds(j,:) = global_inds;
        
        for r=1:size(global_inds,2)
            F_global_external(global_inds(r))= F_global_external(global_inds(r))+ Fgamma(r);
        end
    end
        count = 1;
        del_u = 50;
         while(norm(del_u) > tolerance)     
            K_global = zeros(2*nn,2*nn);
            F_global_internal = zeros(2*nn,1);
            for m = 1:nele
                ele_dof = ele_number(m,:);
                ele_coord = n_coord(ele_number(m,:),:);
                x1 = ele_coord(:,1);
                y1 = ele_coord(:,2);
                f_internal = zeros(8,1);
                K_local = zeros(8,8);
                u_ele = delta_defl(inds(m,:));
                for cc = 1:Np1
                    for pp = 1:Np1
                        eps = Ip1(cc);
                        eta = Ip1(pp);
                        
                        [dum_N,dum_J,B] = shape(eps, eta, x1, y1);
                        d_total = B*u_ele;
                        F = [d_total(1,1) d_total(3,1) 0; d_total(4,1) d_total(2,1) 0; 0 0 0] + kron_delta;
                        C = F'*F;
                        E = 0.5.*(C - kron_delta);
                        [D,T,Sigma] = DTS(F,E,lamda,mu);
                        f_internal = f_internal + B'*Sigma*det(dum_J)*w1(cc)*w1(pp);
                        K_local = K_local + B'*(D+T)*B*det(dum_J)*w1(cc)*w1(pp);
                    end
                end
                
                global_inds = reshape(dof_inds(ele_dof),1,[]); 
                inds(m,:) = global_inds;
                for r=1:size(global_inds,2)
                    for c=1:size(global_inds,2)           
                        K_global(global_inds(r),global_inds(c))= K_global(global_inds(r),global_inds(c))+ K_local(r,c);
                    end 
                    F_global_internal(global_inds(r))= F_global_internal(global_inds(r))+ f_internal(r);
                end
            end
            
            F_global = F_global_external - F_global_internal;
          
            
            n_bc = 1;
            for bb = 1:nn
                if mod(bb,n_x) == 1
                    g1(1,n_bc) = bb;
                    n_bc = n_bc + 1;
                end
            end
            
            dirichlet = reshape(dof_inds(g1),1,[]);
            u_new = zeros(length(dirichlet),1);
            K_global_modified = K_global;
            F_global_modified = F_global;
            
            for dk = 1:size(dirichlet,2)
                K_global_modified(dirichlet(dk),:) = 0.0;
            end
            
            for w=1:size(dirichlet,2)
                F_global_modified(dirichlet(w))= u_new(w);
            end
            
            for pk = 1:size(dirichlet,2)
                F_global_modified(:) = F_global_modified(:) - K_global_modified(:, dirichlet(pk))*u_new(pk);
            end
            
            for oo = 1:size(dirichlet,2)
                K_global_modified(:,dirichlet(oo)) = 0.0;
                K_global_modified(dirichlet(oo),dirichlet(oo)) = 1.0;
            end
          
            del_u = K_global_modified\F_global_modified; 
            delta_defl = delta_defl + del_u;
            count = count +1;
                
         end
   pre = ['Load increment number ',num2str(i),' has converged with ',num2str(count),' iterations'];
   disp(pre)
   center_dof = dof_inds(center_node);
   center_dof_y = center_dof(2);
   def_center_y(i,1) = delta_defl(center_dof_y);
   center_load_P(i,1) = Pnew*l_y;
   load_counter(i,1) = count;
        
end

defl_x_RI = delta_defl(1:2:end);
defl_y_RI = delta_defl(2:2:end);
dof_half_y = n_x*(ele_y/2)+1:1:n_x*(ele_y/2 +1);
defl_half_y = defl_y_RI(dof_half_y);
coord_half_x = n_coord(dof_half_y);


% Stress Calculation and Plotting
i = 1;
for j = 1:ele_y
    ele_left = (j-1/2)*ele_x;
%     coordinate of that element
    coord_ele_left = n_coord(ele_number(ele_left,:),:);
    coord_left_x = coord_ele_left(:,1);
    coord_left_y = coord_ele_left(:,2);
    disp_ele_left = delta_defl(inds(ele_left,:));
    eps = 0.99; eta1 = -0.9; eta2 = 0.9;
    stress_1 = calstress(eps,eta1, coord_left_x, coord_left_y, disp_ele_left, Csmat, kron_delta, lamda, mu);
%     corresponding y-coordinate for plotting
    [Nnew, Jun,Bun] = shape(eps,eta1, coord_left_x, coord_left_y);
    Nmod1 = [Nnew(2,2), Nnew(2,4), Nnew(2,6), Nnew(2,8)];
    ymod1 = Nmod1*coord_left_y;
    
    stress_2 = calstress(eps,eta2, coord_left_x, coord_left_y, disp_ele_left, Csmat, kron_delta, lamda, mu);
%     corresponding y-coordinate for plotting
    [Nnew2, Jun, Bun] = shape(eps, eta2, coord_left_x, coord_left_y);
    Nmod2 = [Nnew2(2,2), Nnew2(2,4), Nnew2(2,6), Nnew2(2,8)];
    ymod2 = Nmod2*coord_left_y;
    
%     Stresses
    sigma_xx(i,1) = stress_1(1);
    sigma_xy(i,1) = stress_1(4);
    corr_coord_half_x(i,1) = ymod1;
    i = i+1;
    sigma_xx(i,1) = stress_2(1);
    sigma_xy(i,1) = stress_2(4);
    corr_coord_half_x(i,1) = ymod2;
    i = i+1;
end

figure(1)
plot(def_center_y,center_load_P)
xlabel('Deflection in Y-direction at D/2');
ylabel('Load P at center (D/2) in lbs');
title('Load vs Deflection with Reduced Integration')

figure(2)
plot(coord_half_x, defl_half_y)
xlabel('x-cordinate')
ylabel('Deflection')
title('Deflection along X-direction @ y = 0 with RI')

figure(3)
plot(corr_coord_half_x, sigma_xx)
xlabel('y-coordinate')
ylabel('Sigma XX in kips')
title('Sigma XX across x = L/2 with RI')

figure(4)
plot(corr_coord_half_x, sigma_xy)
xlabel('y-coordinate')
ylabel('Sigma XY in kips')
title('Sigma XY across x = L/2 with RI')


figure(5)
plotbeam(defl_x_RI,defl_y_RI,ele_x,ele_y)

    
    
    
    
    
    
    
    
    
    










                        
                
 
    
    

                
                
                
                
                
        
    
