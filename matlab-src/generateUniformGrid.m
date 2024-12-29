function [xgrid,xc,dx] = generateUniformGrid(Nx,Lx)
%Generates a uniformly spaced grid
%dx - cell sizes
%xgrid - cell node locations
%xc - cell centre locations
xx = (0:Nx)/Nx;
xgrid(1:Nx+1) = xx(1:Nx+1)*Lx;
dx = xgrid(2) - xgrid(1);
xc(1:Nx) = 0.5*(xgrid(1:Nx)+xgrid(2:Nx+1));
end

