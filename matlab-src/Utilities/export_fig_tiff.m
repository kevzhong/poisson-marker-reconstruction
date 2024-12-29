function export_fig_tiff(fh,filename, append)
if nargin == 1
    append = 0;
end

if endsWith(filename, '.tiff')
    filename(end-4:end) = [];
end

% if ~append
%     fprintf('exporting tiff... %s\n', filename);
%     export_fig(filename,'-tiff','-m2','-nocrop');
%     fprintf('\tdone.\n');
% elseif append == 1
%     fprintf('appending to tiff... %s\n', filename);
    export_fig(fh,filename,'-tiff','-m2','-nocrop','-append');
    
% end

end