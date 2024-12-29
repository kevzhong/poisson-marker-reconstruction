function setLabels(legendFont,legendLoc,xlab,ylab,xlims,ylims)
%Formats a current plot


if ~isempty(legendFont) && ~isempty(legendLoc)
    legend({},'Interpreter','Latex','Fontsize',legendFont,'location',legendLoc,'box','off')
end
xlabel(xlab)
ylabel(ylab)
%set(gca,'TickLabelInterpreter','latex')

if ~isempty(xlims)
    xlim(xlims);
end

if ~isempty(ylims)
    ylim(ylims);
end
end

