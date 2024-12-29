function print_fig(fHandle,fpath,fType)

set(fHandle,'paperpositionmode','auto')

if isequal(fType,'-png')
    print(fHandle,fpath,'-dpng','-painters');
elseif isequal(fType,'-pdf')
    print(fHandle,fpath,'-depsc','-painters');
    cmd = strcat('epstopdf',{' '},fpath,'.eps',{' '},fpath,'.pdf');
    cmd = cmd{1};
    system(cmd);
    %cmd = strcat('rm',{' '},fpath,'.eps');
    %cmd = cmd{1};
    %system(cmd);
else
    error('Unrecognised file type. Must be either ''-png'' or ''-pdf''. ')
end