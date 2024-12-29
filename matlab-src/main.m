clean
init;


f = figure('color','white','Units','centimeters');

fwidth = 32;
gap = [0.1 0.11];
marg_h = [0.11 0.11];
marg_w = [0.11 0.11];
[ax,~] = tight_subplot(1,2,gap,marg_h,marg_w);
arrayfun(@(x) hold(x, 'on'), ax);
arrayfun(@(x) box(x, 'on'), ax);
arrayfun(@(x) set(x,'XTickLabelMode','auto'),ax);
arrayfun(@(x) set(x,'YTickLabelMode','auto'),ax);
arrayfun(@(x) set(x,'FontSize',16),ax);
arrayfun(@(x) set(x,'Layer','top'),ax);
arrayfun(@(x) set(x,'TickLength',x.TickLength*3),ax);
arrayfun(@(x) daspect(x,[1 1 1]),ax);

f.Position(2) = 10;
f.Position(3) = fwidth;
f.Position(4) = fwidth/2;

% Given a Lagrangian representation of a geometry
% specified by panels of size ds, oriented normals nhat
% Construct a [0,1] indicator function field, where the geometry is
% located by the 0.5 contour

% Eulerian grid points
Nx = 128; 
Ny = 128;
% Lagrangian node count
NL = 188; 

% %Eulerian grid points
% Nx = 64; 
% Ny = 64;
% %Lagrangian node count
% NL = 94; 



Lx = 5;
Ly = 5;

[~,xm,dx] = generateUniformGrid(Nx,Lx);
[~,ym,dy] = generateUniformGrid(Ny,Ly);

xc = xm - dx/2;
yc = ym - dy/2;

%%%%%%%%%%%% GENERATE LAGRANGIAN MESH %%%%%%%%%%%%%%%%%%%%%%

% % Unit circle
% R = 1;
% theta = linspace(0,2*pi,NL+1);
% xv = R*cos(theta);
% yv = R*sin(theta);
% DI = -1.0; % Interface jump value


% Star-shape
R1 = 1;
R2 = 0.2*R1;
theta = linspace(0,2*pi,NL+1);
n=5;
r = R1 + R2 * sin(n*theta + pi/2);
xv = r.*sin(theta);
yv = r.*cos(theta);
DI = 1.0; % Interface jump value


xv = xv + Lx/2;
yv = yv + Ly/2;

xl = 0.5* (xv(1:end-1) + xv(2:end) );
yl = 0.5* (yv(1:end-1) + yv(2:end) );
dyl =  yv(2:end) - yv(1:end-1) ;
dxl = xv(2:end)  - xv(1:end-1) ;
ds = sqrt(    dxl.^2   +   dyl.^2   );
nhat = zeros(2,numel(ds));

nhat(1,:) = dyl;
nhat(2,:) = -dxl;

nhat = nhat ./ vecnorm(nhat,2,1);



plot(ax(1),xv,yv,'k-')
quiver(ax(1),xl,yl,-nhat(1,:),-nhat(2,:),0.5,'k-')
plot(ax(1),xl,yl,'bs','MarkerFaceColor','w','MarkerSize',4)
xlim(ax(1),[0 Lx]) ; ylim(ax(1),[0 Ly]);
title(ax(1),'Lagrangian front','Fontweight','normal')

%error('asdf')
ds_on_dx = mean(ds)/dx


% Gradient of indicator function
Gx = zeros(Nx,Ny);
Gy = zeros(Nx,Ny);

%dfunc = @(r) romaDelta(r);
%dfunc = @(r) pmDelta(r);

%dfunc = @(r) peskinDelta(r);
dfunc = @(r) brDelta(r);

nsup = 2;

% Begin Lagrangian calculation for gradients
for n = 1:NL

    %%%%%%%%%%%%%%%%    X-gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Staggered Eulerian indices of nearest neighbour containing node
    i = round(xl(n) / dx ) + 1; 
    j = floor(yl(n) / dy ) + 1;

    %j = round(yl(n) / dy ) + 1;

        for ii = i-nsup:i+nsup
        rxc = ( xl(n) - xc(ii) ) / dx;
        drxc = dfunc(rxc);
        for jj = j-nsup:j+nsup
            rym = ( yl(n) - ym(jj) ) / dy;
            drym = dfunc(rym);
            w_imh_j = drxc * drym; %x-stagger
            Gx(ii,jj) = Gx(ii,jj) + DI * nhat(1,n)  * w_imh_j * ds(n) / (dx*dy);
        end
        end

        %%%%%%%%%%%%%%%%    Y-gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        i = floor(xl(n) / dx ) + 1;
        j = round(yl(n) / dy ) + 1;


    % 3x3 interpolation domain
    for ii = i-nsup:i+nsup

        rxm = ( xl(n) - xm(ii) ) / dx;
        drxm = dfunc(rxm);

        for jj = j-nsup:j+nsup

            ryc = ( yl(n) - yc(jj) ) / dy;

            dryc = dfunc(ryc);

            w_i_jmh = drxm * dryc; %y-stagger

            Gy(ii,jj) = Gy(ii,jj) + DI * nhat(2,n)  * w_i_jmh * ds(n) / (dx*dy);

        end
    end

end



%%%%%%%%%%   Begin solution to Poisson equation %%%%%%%%%%%%%%%%%%%%%%%
rhs = zeros(Nx,Ny);

for i = 1:Nx-1
    for j = 1:Ny-1

        rhs(i,j) = ( Gx(i+1,j) - Gx(i,j) ) / dx   + ...
                   ( Gy(i,j+1) - Gy(i,j) ) / dy ;
    end
end



% %%%%%%%%%%%%%%%%%%%%%%   FAST-POISSON-SOLVE %%%%%%%%%%%%%%%%%%%%%%%%%%

TYPE = 'FFT';
%TYPE = 'FFT';

%Wavenumber vectors
kx = 0:(Nx-1);
ky = 0:(Ny-1);

switch TYPE
    case 'DCT'
        % % DCT Modified wavenumbers
        lmb_x_on_dx2 = 2* (cos( pi*kx./Nx ) - 1);
        lmb_y_on_dy2 = 2* (cos( pi*ky./Ny ) - 1);
    case 'FFT'
        % FFT Modified wavenumbers
        lmb_x_on_dx2 = 2* (cos( 2*pi*kx./Nx ) - 1);
        lmb_y_on_dy2 = 2* (cos( 2*pi*ky./Ny ) - 1);
    otherwise
        error('invalid transform type specified, typo?')
end

lmb_x_on_dx2 = lmb_x_on_dx2 ./ dx^2;
lmb_y_on_dy2 = lmb_y_on_dy2 ./ dy^2;


switch TYPE
    case 'DCT'
        rhs_hat = dct2(rhs) ; %/ sqrt(2 * Nx * Ny); % Type-II implicit
    case 'FFT'
        rhs_hat = fft2(rhs) / (Nx * Ny);
end


vhat = rhs_hat(Nx,Ny);


for k = 1:Nx
    for m = 1:Ny
        vhat(k,m) = rhs_hat(k,m) / ( lmb_x_on_dx2(k) + lmb_y_on_dy2(m)  );
    end
end

% Arbitrary
vhat(1,1) = 0.0;

switch TYPE
    case 'DCT'
        vof = idct2(vhat);
    case 'FFT'
        vof = real( ifft2(vhat) ) ;
end


% Rescale [min,max] range to [0,1]
minvof = min(min(vof));
maxvof = max(max(vof));

vof = (vof - minvof ) ./ (maxvof - minvof) ;

% %%%%%%%%%%%%%%%%%%%%%%   GAUSS-SEIDEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% tolerance = 1e-6;
% N_k       = 100000;
% 
% 
% r = zeros(Nx,Ny);
% vof = zeros(Nx,Ny);
% 
% % Compute initial residual
% for i=2:Nx-1
%     for j=2:Ny-1
%         r(i,j) = -rhs(i,j) + ...
%                  ( vof(i-1,j) - 2.0 * vof(i,j) + vof(i+1,j) ) / dx^2 + ...
%                  ( vof(i,j-1) - 2.0 * vof(i,j) + vof(i,j+1) ) / dy^2 ;
%     end
% end
% r_norm	= sqrt(sum(sum(r.^2)));
% 
% k = 0;
% % Gauss-Seidel iterative loop
% while r_norm>tolerance && k<N_k
% 
%     % Gauss-Seidel iteration
%     for i=2:Nx-1
%         for j=2:Ny-1
%             vof(i,j) = (dx^2*-rhs(i,j) + vof(i-1,j) + vof(i+1,j) + vof(i,j-1) + vof(i,j+1))/4;
%         end
%     end
% 
%     % Compute residual
%     for i=2:Nx-1
%         for j=2:Ny-1
%             r(i,j) = dx^2*-rhs(i,j) + vof(i-1,j) + vof(i+1,j) + vof(i,j-1) + vof(i,j+1) - 4*vof(i,j);
%         end
%     end
% 	r_norm = sqrt(sum(sum(r.^2)));
%     k      = k + 1
% 
% 
% 
% end
% %%%%%%%%%%%%%%%%%%% END GAUSS-SEIDEL %%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%vof = abs(vof);

mean_VOF = mean(vof,'all')
%area = pi * R^2

%surf(ax(2),xc',ym,vof')

cmap = brewermap(256,'RdBu');

%p = pcolor(ax(2),xc',ym,vof');
%p.EdgeColor = 'none';
hold(ax(2),'on');
imagesc(ax(2),xc',ym,vof')
xlim(ax(2),[0 Lx]) ; ylim(ax(2),[0 Ly]);
%contour(ax(2),xc',ym,vof',[0.5 0.5],'b-','Linewidth',2)
daspect(ax(2), [1 1 1] );
axPos = ax(2).Position;
colorbar(ax(2));
clim(ax(2),[0 1]);
colormap(ax(2),cmap)
ax(2).Position = axPos;

title(ax(2),'Reconstructed $\phi$','Interpreter','latex','Fontweight','normal')

