function [h_main, h_inset,main_legend]=inset(main_handle, inset_handle,inset_width,inset_height,location)
% The function plotting figure inside figure (main and inset) from 2 existing figures.
% inset_size is the fraction of inset-figure size, default value is 0.35
% The outputs are the axes-handles of both.
%
% An examle can found in the file: inset_example.m
%
% Moshe Lindner, August 2010 (C).
% if nargin==2
%     inset_size=0.35;
% end
inset_width=inset_width*.7;
inset_height=inset_height*.7;

figure
new_fig=gcf;
main_fig = findobj(main_handle,'Type','axes');
main_legend = findobj(main_handle,'Type','Legend');
h_main = copyobj(main_fig,new_fig);
set(h_main,'Position',get(main_fig,'Position'))
inset_fig = findobj(inset_handle,'Type','axes');
h_inset = copyobj(inset_fig,new_fig);
ax=get(main_fig,'Position');
switch location
    %ax(1) left
    %ax(2) bottom
    %ax(3) width
    %ax(4) height
    case 'northeast'
        set(h_inset,'Position', [.7*ax(1)+ax(3)-inset_width .9*ax(2)+ax(4)-inset_height inset_width inset_height])
    case 'northwest'
        set(h_inset,'Position', [ax(1) + 0.55*ax(1) .9*ax(2)+ax(4)-inset_height inset_width inset_height])
    case 'southwest'
        set(h_inset,'Position', [ax(1) + 0.55*ax(1) 1.4*ax(2)  inset_width inset_height])
    case 'southeast'
        set(h_inset,'Position', [.7*ax(1)+ax(3)-inset_width 1.4*ax(2)  inset_width inset_height])
    otherwise
        error('Invalid inset location.')
end
yticks(0:1:6)
%set(h_main,'Legend',main_legend);
end