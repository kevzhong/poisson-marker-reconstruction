clean
init;

DI = -1.0;


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

view(ax(1),-45,35.264)

f.Position(2) = 10;
f.Position(3) = fwidth;
f.Position(4) = fwidth/2;

% Given a Lagrangian representation of a geometry
% specified by panels of size ds, oriented normals nhat
% Construct a [0,1] indicator function field, where the geometry is
% located by the 0.5 contour

% Nx = 64; 
% Ny = 64;
% Nz = 64;
% fname = 'gts/unitSphere_N3.gts';

% % % Eulerian grid points
% Nx = 128; 
% Ny = 128;
% Nz = 128;
% fname = 'gts/unitSphere_N4.gts';
% Lx = 4;
% Ly = 4;
% Lz = 4;

% % Eulerian grid points
Nx = 128; 
Ny = 128;
Nz = 128;
fname = 'gts/ABC_ellipsoid_N32.gts';
Lx = 8;
Ly = 8;
Lz = 8;


% % % Eulerian grid points
% Nx = 128; 
% Ny = 128;
% Nz = 128;
% fname = 'gts/alpha_6k.gts';
% Lx = 1.5;
% Ly = 1.5;
% Lz = 1.5;


% % % Eulerian grid points
% Nx = 320; 
% Ny = 320;
% Nz = 320;
% fname = 'gts/bunny.gts';
% Lx = 0.2;
% Ly = 0.2;
% Lz = 0.2;

[~,xm,dx] = generateUniformGrid(Nx,Lx);
[~,ym,dy] = generateUniformGrid(Ny,Ly);
[~,zm,dz] = generateUniformGrid(Nz,Lz);

xc = xm - dx/2;
yc = ym - dy/2;
zc = zm - dz/2;


% Read triangulated geometry
[xyzv,xyz,verts_of_face,nhat,Atri] = read_gts(fname);


% Re-Centre
xyzv(:,1) = xyzv(:,1) + Lx/2;
xyzv(:,2) = xyzv(:,2) + Ly/2;
xyzv(:,3) = xyzv(:,3) + Lz/2;

xyz(:,1) = xyz(:,1) + Lx/2;
xyz(:,2) = xyz(:,2) + Ly/2;
xyz(:,3) = xyz(:,3) + Lz/2;


% % Bunny re-centre
% xyzv(:,1) = xyzv(:,1) + 0.12;
% xyzv(:,2) = xyzv(:,2) ;
% xyzv(:,3) = xyzv(:,3) + Lz/2;
% 
% xyz(:,1) = xyz(:,1) + 0.12;
% xyz(:,2) = xyz(:,2) ;
% xyz(:,3) = xyz(:,3) + Lz/2;

NL = size(nhat,1);

Atri_on_dx2 = mean(Atri) / dx^2



% Bouding box
xmin = min( xyzv(:,1) );
xmax = max( xyzv(:,1) );
DLX = xmax - xmin
ymin = min( xyzv(:,2) );
ymax = max( xyzv(:,2) );
DLY = ymax - ymin
zmin = min( xyzv(:,3) );
zmax = max( xyzv(:,3) );
DLZ = zmax - zmin
%error('asdf')

p1 = patch(ax(1),'faces',verts_of_face,'vertices',xyzv,'FaceColor',[0.8 0.8 0.8],'Linewidth',0.5); %camlight
%p1 = patch(ax(1),'faces',verts_of_face,'vertices',xyzv,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); camlight(ax(1));
xlabel(ax(1),'x') ; ylabel(ax(1),'y') ; zlabel(ax(1),'z');

%error('asdf')

dfunc = @(r) brDelta(r);
%dfunc = @(r) romaDelta(r);
%dfunc = @(r) pmDelta(r);
%dfunc = @(r) peskinDelta(r);

nsup = 2;

% Gradient of indicator function
Gx = zeros(Nx,Ny,Nz);
Gy = zeros(Nx,Ny,Nz);
Gz = zeros(Nx,Ny,Nz);

for n = 1:NL

    xl = xyz(n,1) ; yl = xyz(n,2) ; zl = xyz(n,3);


    %%%%%%%%%%%%%%%%    X-gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Staggered Eulerian indices of nearest neighbour containing node
    i = round(xl / dx ) + 1;
    j = floor(yl / dy ) + 1;
    k = floor(zl / dz ) + 1;

    for ii = i-nsup:i+nsup
        rxc = ( xl - xc(ii) ) / dx;
        drxc = dfunc(rxc);
        for jj = j-nsup:j+nsup
            rym = ( yl - ym(jj) ) / dy;
            drym = dfunc(rym);
            for kk = k-nsup:k+nsup
                rzm = ( zl - zm(kk) ) / dz;
                drzm = dfunc(rzm);
                weight = drxc * drym * drzm; %x-stagger
                Gx(ii,jj,kk) = Gx(ii,jj,kk) + DI * nhat(n,1)  * weight * Atri(n) / (dx*dy*dz);
            end
        end
    end


    %%%%%%%%%%%%%%%%    Y-gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Staggered Eulerian indices of nearest neighbour containing node
    i = floor(xl / dx ) + 1;
    j = round(yl / dy ) + 1;
    k = floor(zl / dz ) + 1;

    for ii = i-nsup:i+nsup
        rxm = ( xl - xm(ii) ) / dx;
        drxm = dfunc(rxm);
        for jj = j-nsup:j+nsup
            ryc = ( yl - yc(jj) ) / dy;
            dryc = dfunc(ryc);
            for kk = k-nsup:k+nsup
                rzm = ( zl - zm(kk) ) / dz;
                drzm = dfunc(rzm);
                weight = drxm * dryc * drzm; %x-stagger
                Gy(ii,jj,kk) = Gy(ii,jj,kk) + DI * nhat(n,2)  * weight * Atri(n) / (dx*dy*dz);
            end
        end
    end


    %%%%%%%%%%%%%%%%    Z-gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Staggered Eulerian indices of nearest neighbour containing node
    i = floor(xl / dx ) + 1;
    j = floor(yl / dy ) + 1;
    k = round(zl / dz ) + 1;

    for ii = i-nsup:i+nsup
        rxm = ( xl - xm(ii) ) / dx;
        drxm = dfunc(rxm);
        for jj = j-nsup:j+nsup
            rym = ( yl - ym(jj) ) / dy;
            drym = dfunc(rym);
            for kk = k-nsup:k+nsup
                rzc = ( zl - zc(kk) ) / dz;
                drzc = dfunc(rzc);
                weight = drxm * drym * drzc; %x-stagger
                Gz(ii,jj,kk) = Gz(ii,jj,kk) + DI * nhat(n,3)  * weight * Atri(n) / (dx*dy*dz);
            end
        end
    end


end

disp('Finished Lagrangian loop')

%imagesc(ax(2),squeeze(Gz(:,:,Nz/2) )' )
% Build RHS to Poisson equation

%%%%%%%%%%   Begin solution to Poisson equation %%%%%%%%%%%%%%%%%%%%%%%
rhs = zeros(Nx,Ny,Nz);
for i = 1:Nx-1
    for j = 1:Ny-1
        for k = 1:Nz-1
            rhs(i,j,k) = ( Gx(i+1,j,k) - Gx(i,j,k) ) / dx   + ...
                       ( Gy(i,j+1,k) - Gy(i,j,k) ) / dy   + ...
                       ( Gz(i,j,k+1) - Gz(i,j,k) ) / dz   ;
        end
    end
end

% %%%%%%%%%%%%%%%%%%%%%%   FAST-POISSON-SOLVE %%%%%%%%%%%%%%%%%%%%%%%%%%

TYPE = 'FFT';

%Wavenumber vectors
kx = 0:(Nx-1);
ky = 0:(Ny-1);
kz = 0:(Nz-1);

switch TYPE
    case 'DCT'
        % % DCT Modified wavenumbers
        lmb_x_on_dx2 = 2* (cos( pi*kx./Nx ) - 1);
        lmb_y_on_dy2 = 2* (cos( pi*ky./Ny ) - 1);
        lmb_z_on_dz2 = 2* (cos( pi*kz./Nz ) - 1);

    case 'FFT'
        % FFT Modified wavenumbers
        lmb_x_on_dx2 = 2* (cos( 2*pi*kx./Nx ) - 1);
        lmb_y_on_dy2 = 2* (cos( 2*pi*ky./Ny ) - 1);
        lmb_z_on_dz2 = 2* (cos( 2*pi*kz./Nz ) - 1);

    otherwise
        error('invalid transform type specified, typo?')
end

lmb_x_on_dx2 = lmb_x_on_dx2 ./ dx^2;
lmb_y_on_dy2 = lmb_y_on_dy2 ./ dy^2;
lmb_z_on_dz2 = lmb_z_on_dz2 ./ dz^2;

disp('Solving Poisson...')

switch TYPE
    case 'DCT'
        rhs_hat = dct(rhs,[],1);
        rhs_hat = dct(rhs_hat,[],2);
        rhs_hat = dct(rhs_hat,[],3);

        %rhs_hat = dct2(rhs) ; %/ sqrt(2 * Nx * Ny); % Type-II implicit
    case 'FFT'
        %rhs_hat = fft2(rhs) / (Nx * Ny);
        rhs_hat = fftn(rhs);
end


vhat = rhs_hat;


for k = 1:Nx
    for m = 1:Ny
        for n = 1:Nz
            vhat(k,m,n) = rhs_hat(k,m,n) / ( lmb_x_on_dx2(k) + lmb_y_on_dy2(m) + lmb_z_on_dz2(n) );
        end
    end
end

% Arbitrary
vhat(1,1,1) = 0.0;

switch TYPE
    case 'DCT'
        %vof = idct2(vhat);
        vof = idct(vhat,[],1);
        vof = idct(vof,[],2);
        vof = idct(vof,[],3);

    case 'FFT'
        vof = real( ifftn(vhat) ) ;
end


% Rescale [min,max] range to [0,1]
minvof = min(vof,[],'all');
maxvof = max(vof,[],'all');

vof = (vof - minvof ) ./ (maxvof - minvof) ;

imagesc(ax(2),xm',ym,squeeze(vof(:,Ny/2,:) )' )

cmap = brewermap(256,'RdBu');
daspect(ax(2), [1 1 1] );
axPos = ax(2).Position;
colorbar(ax(2));
clim(ax(2),[0 1]);
colormap(ax(2),cmap)
ax(2).Position = axPos;

title(ax(2),'Reconstructed $\phi$','Interpreter','latex','Fontweight','normal')

Mat2VTK('vof.vtk',vof,'binary',[dx dy dz])

%figure(10)
%plot(squeeze(vof(:,Ny/2,Nz/2) ),'k-o')
