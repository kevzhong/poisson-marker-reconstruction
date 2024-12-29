function formatPlot(ax,legendLoc,xlab,ylab,xlims,ylims,latex_flag)
%Formats a current plot
% lab_size = 11;
% ticklabSize = 11;
% legendFont = 11;

% lab_size = 24;
% ticklabSize = 24;
% legendFont = 24;

if nargin == 6
    latex_flag = true; %Assume using latex if flag not specified
end
%gcf;
if isempty(ax)
    ax = gca;
end
%ax.FontSize = ticklabSize;

if ~isempty(legendLoc)
    legend(ax,{},'location',legendLoc,'box','off')
    if latex_flag
            ax.Legend.Interpreter = 'latex';
    end
end
xlabel(ax,xlab)%,'Interpreter','Latex')
ylabel(ax,ylab)%,'Interpreter','Latex')
%set(gca,'TickLabelInterpreter','latex')

if ~isempty(xlims)
    xlim(ax,xlims);
end

if ~isempty(ylims)
    ylim(ax,ylims);
end

if latex_flag
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
end

