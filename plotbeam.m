function plotbeam(ux,uy,ele_x,ele_y)
l_x = 50; 
l_y = 5; 
a_x = l_x/ele_x;
a_y = l_y/ele_y; 
mf = 1;
mu_x = reshape(ux,ele_x + 1,ele_y +1)';

mu_y = reshape(uy,ele_x+1,ele_y+1)';

[grid_x,grid_y] = meshgrid(0:a_x:l_x,0:a_y:l_y);

% plot the solution on the deformed configuration

surf(grid_x+ (mf)*mu_x, grid_y + (mf*mu_y) , (mf*mu_y), 'FaceColor', 'interp')

colormap jet; 
view(2); 
end