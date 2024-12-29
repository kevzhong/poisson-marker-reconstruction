function rdelta = romaDelta(r)
% Compute the Roma delta function


if (r <= 1.5) && (r >= 0.5)
    rdelta = 5.0 - 3.0*abs(r) - sqrt(-3.0*(1.0 - abs(r))^2 + 1.0);
    rdelta = rdelta / 6.0;
elseif (r <= 0.5) && (r >= 0.0)
    rdelta = 1.0 + sqrt(-3.0*r^2 + 1.0);
    rdelta = rdelta / 3.0;
else
    rdelta = 0.0;
end