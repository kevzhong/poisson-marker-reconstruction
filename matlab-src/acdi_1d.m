clean

% Simple constant 1D advection problem

N = 128; %No. of spatial grid points
nt = 100; %No. of timesteps

tframe = 1;

Uc = 5.0; % Advection velocity

small = 1e-15;


%dt = 1e-5;
dt = 1e-4;

%%%% Plotting axes limits %%%%%%%%%
limflag = true;


%%%%%%%%%%%%%  Generate grid %%%%%%%%%%%%%%%%%%%%%%%
L = 1;
[~,x,dx] = generateUniformGrid(N,L);



%%%%%%%%%%%%%  Create initial conditions %%%%%%%%%%%
R0 = 0.2;
%epsilon = 0.55*dx;
epsilon = dx;


phi_init = zeros(1,N);
small2 = small;

for k = 1:N
    rr = sqrt ( ( x(k) - L/2 )^2 );
    tmp = 0.5 * (1 - tanh( (rr - R0 ) /(2*epsilon) )  );
    tmp = max(0.0 + small2, tmp);
    tmp = min(1.0 - small2, tmp);
    phi_init(k) = tmp;
end


%RK3 coefficients
gamma = [8/15,5/12,3/4];
zeta = [0,-17/60,-5/12];
alp = gamma + zeta;


% fwidth = 24;
% fheight = fwidth * 1.4;

fwidth = 48;
fheight = 32;


f = figure('color','white','Units','centimeters');
f.Position(3) = fwidth;
f.Position(4) = fheight;
gap = [0.05 0.1];
marg_h = [0.1 0.03];
marg_w = [0.2 0.2];
[ha,~] = tight_subplot(3,1,gap,marg_h,marg_w);

arrayfun(@(x) hold(x,'on'), ha);
arrayfun(@(x) box(x, 'on'), ha);
arrayfun(@(x) set(x,'XTickLabelMode','auto'),ha);
arrayfun(@(x) set(x,'YTickLabelMode','auto'),ha);
arrayfun(@(x) set(x,'FontSize',20),ha);
arrayfun(@(x) set(x,'Layer','top'),ha);
arrayfun(@(x) set(x,'TickLength',x.TickLength*3),ha);
%f.Visible = 'off';



phi = phi_init;
u = ones(size(phi)) .* Uc;

%non-linear terms at previous RK3 substep
%reset at RK3 beginning
rhs_m1 = zeros(size(phi));
rhs = zeros(size(phi));


psi = zeros(size(phi));

nhat = zeros(size(phi));

kmv = (1:N) - 1; kmv(1) = N;
kpv = (1:N) + 1; kpv(N) = 1;

gamACDI = abs( Uc ) ;



%psi = sgnd(eps,phi);

cnt = 0;


for i = 2:nt

    %for nss = 1:3
    %    gamdt = gamma(nss) * dt;
    %    zetdt = zeta(nss)  * dt;


    % Compute signed distance function
    for k = 1:N

        small2 = 1e-15;


        tmp = max(phi(k),0 + small2);
        tmp = min(tmp, 1 - small2);
        psi(k) = epsilon * log( (tmp + small2) / ( 1 - tmp + small2 )  );

        %psi(k) = epsilon * log( (phi(k) + small2) / ( 1 - phi(k) + small2 )  );
    end

    % Compute normals
    for k = 1:N
        km = kmv(k) ; kp = kpv(k);
        dpsi = 0.5 * ( psi(kp) - psi(km) ) / dx;


        if abs(dpsi) > small
            inv_magN = 1.0 / ( sqrt(dpsi^2) + small );
            nhat(k) = dpsi * inv_magN;
        else
            nhat(k) = 0;
        end

    end

    % % Compute RHS fluxes: ph and mh collected separately
    for k = 1:N
        km = kmv(k) ; kp = kpv(k);

        %Advection
        rhs(k) =  -(u(kp)*0.5*(phi(kp)+phi(k) ) - u(k)*0.5*(phi(k)+phi(km))) / dx ;

        %Diffusive
        rhs(k) = rhs(k) + gamACDI*epsilon*(phi(kp)-2*phi(k)+phi(km)) / (dx^2);


        % Anti-diffusion
        psi_ph = 0.5*(psi(k ) + psi(kp) );
        psi_mh = 0.5*(psi(km) + psi(k ) );
        nhat_ph = 0.5*(nhat(k ) + nhat(kp) );
        nhat_mh = 0.5*(nhat(km) + nhat(k ) );

        f_kp = -0.25*(1.0 - tanh(0.5*psi_ph/epsilon )^2 * nhat_ph );
        f_km = -0.25*(1.0 - tanh(0.5*psi_mh/epsilon )^2 * nhat_mh );


        rhs(k) = rhs(k) + gamACDI * 0.5 * (f_kp - f_km) / dx;


    end



    phi = phi + dt * rhs;

    phimin = min(phi);
    phimax = max(phi);


    if ( mod(i,tframe) == 0)
        cla(ha(1));
        str = sprintf('$tU/R_0 = %.8f$',(i-1)*dt*Uc/R0);
        title(ha(1),str,'Interpreter','latex','FontWeight','bold');
        plot(ha(1),x,phi,'k-','Linewidth',2)
        ylabel(ha(1),'$\phi(x,t)$','Interpreter','latex')
        ha(1).XTickLabel = [];


        cla(ha(2));
        title(ha(2),str,'Interpreter','latex','FontWeight','bold');
        plot(ha(2),x,psi/epsilon,'b-','Linewidth',2)
        yline(ha(2),0,'k--');

        ha(2).XTickLabel = [];
        ylabel(ha(2),'$\psi(x,t)/\epsilon$','Interpreter','latex')
        %if limflag
        %    ylim(ha(2),lims_c);
        %end
        drawnow

        cla(ha(3));
        plot(ha(3),x,nhat,'r-','Linewidth',2)
        %xline(ha(3),Dcent - R0,'k--');
        %xline(ha(3),Dcent + R0,'k--');

        xlabel(ha(3),'x')
        ylabel(ha(3),'$\hat{n}$','Interpreter','latex')
        ylim(ha(3),[-2 2])

        %yline(ha(3),-activity,'k--');
        %yline(ha(3),activity,'k--');

        %if limflag
        %    ylim(ha(3),lims_dcdn);
        %end
        drawnow

        cnt = cnt + 1;
        str = sprintf('%08d',cnt); %Off-by-one index
        export_fig(gcf,['~/surfdrive/temp/acdi_1d/frame',str,'.png'],'-m3','-nocrop');
    end



    disp(i)

end


%plot(x,phi,'r-','Linewidth',2)
%plot(x,psi,'b-','Linewidth',2)
%plot(x,nhat,'r-','Linewidth',2)

