function dz_ghost = get_dz_ghost(dz)

dz_ghost = [dz(3:-1:1); dz; dz(end:-1:end-2)];


end