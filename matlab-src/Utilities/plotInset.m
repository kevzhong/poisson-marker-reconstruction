function plotInset(fighandle,x,y,plotType,location,label_x,label_y,lc)
%Plots (x,y) as an inset in the target figure handle

% normalized units are used for position
% parent - handle to main figure
% position - [x,y,width,height] of inset

figwidth = 0.3;
figheight = 0.3;
%lc = 'k';
lw = 2;

switch location %Left, Bottom, Width, Height
    case 'northeast'
        axes('parent',fighandle,'position',[0.9 - figwidth, 0.9-figheight, figwidth, figheight]);
    case 'northwest'
        axes('parent',fighandle,'position',[1-figwidth-0.5, 0.9 - figheight, figwidth, figheight]);
    case 'southwest'
        axes('parent',fighandle,'position',[1-figwidth-0.5, 1-figheight-0.4, figwidth, figheight]);
    case 'southeast'
        axes('parent',fighandle,'position',[0.9 - figwidth, 1-figheight-0.4,  figwidth, figheight]);
end

switch plotType
    case 'semilogx'
        semilogx(x,y,'color',lc,'Linewidth',lw)
        xticks([1,10,100])
    case 'semilogy'
        semilogy(x,y,'color',lc,'Linewidth',lw)
    case 'linear'
        plot(x,y,'color',lc,'Linewidth',lw)
end

ax = gca;
ax.FontSize = 12;
xlabel(label_x,'Fontsize',12,'Interpreter','Latex')
ylabel(label_y,'Fontsize',12,'Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
hold on
%tightfig(gcf);
end

