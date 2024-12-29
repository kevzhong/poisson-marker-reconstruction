function ps_frag(fHandle,fpath,originals,replacements)

assert(length(originals) == length(replacements));


set(fHandle,'paperpositionmode','auto')


print(fHandle,fpath,'-depsc2','-loose','-painters');


%Retrieve file name from fpath
ind = regexp(fpath,'/');
ind = ind(end);

fpath = char(fpath);
fpath_pref = fpath(1:ind);

fname = fpath(ind+1:end);


%Create latex file
fid = fopen(strcat(fpath_pref,fname,'.tex'),'w');
fprintf(fid,'\\documentclass[10pt]{standalone}\n');
fprintf(fid,'\\usepackage{amsmath}\n');
fprintf(fid,'\\usepackage{upgreek}\n');
fprintf(fid,'\\usepackage{psfrag}\n');

fprintf(fid,'\\newcommand\\Rey{\\mbox{\\textit{Re}}}\n');
fprintf(fid,'\\newcommand\\Pran{\\mbox{\\textit{Pr}}}\n');
fprintf(fid,'\\newcommand\\Ray{\\mbox{\\textit{Ra}}}\n');
fprintf(fid,'\\newcommand\\Nus{\\mbox{\\textit{Nu}}}\n');
fprintf(fid,'\\newcommand\\Stan{\\mbox{\\textit{St}}}\n');


fprintf(fid,'\\begin{document}\n');


%PS-frag contents
for i = 1:length(originals)
    fprintf(fid,'\\psfrag{%s}{%s}\n',originals{i},replacements{i});
end
fprintf(fid,'\\includegraphics{%s.eps}\n',fname);
fprintf(fid,'\\end{document}\n');


fclose(fid);


%Execute and convert to pdf
%fid = fopen('')
%cmd = strcat('epstopdf',{' '},fpath,'.eps',{' '},fpath,'.pdf');
%cmd = cmd{1};
%system(cmd);

[~,~] = system(['cd ',fpath_pref, ' && ', 'latex ',fpath_pref,fname,'.tex']);
[~,~] = system(['cd ',fpath_pref, ' && ', 'dvips ',fpath_pref,fname,'.dvi']);
[~,~] = system(['cd ',fpath_pref, ' && ', 'ps2pdf ',fpath_pref,fname,'.ps']);


%File clean-up
f_exts = {'aux','dvi','log','ps','tex'};

for i = 1:length(f_exts)
    delete(strcat(fpath,'.',f_exts{i}));
end


end