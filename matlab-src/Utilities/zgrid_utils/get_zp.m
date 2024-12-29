function zp = get_zp(dz)
zp = (cumsum([0;dz(1:end-1)]) + cumsum(dz) )/2;
