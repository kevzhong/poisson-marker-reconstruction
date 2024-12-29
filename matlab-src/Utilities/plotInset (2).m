function plotInset(fighandle,x,y)
%Plots (x,y) as an inset in the target figure handle

% normalized units are used for position
% parent - handle to main figure
% position - [x,y,width,height] of inset
axes('parent',fighandle,'position',[0.21 0.35 0.3 0.3]);

semilogx(x,y,'k','Linewidth',2)
    ax = gca;
    ax.FontSize = 12;
    xlabel('$z^+$','Fontsize',12,'Interpreter','Latex')
    ylabel('$\Delta U^+$','Fontsize',12,'Interpreter','Latex')
    set(gca,'TickLabelInterpreter','latex')
end

