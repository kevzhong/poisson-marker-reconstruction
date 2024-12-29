function dz_grid = read_dz(dz_dir)

fid = fopen(strcat(dz_dir,'dz_grid.dat'),'r');

dz_grid = fread(fid,'double');

fclose(fid);
end