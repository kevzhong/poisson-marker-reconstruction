function rdelta = brDelta(r)
% Compute the Brackbill & Ruppel (1986) delta function


if ( abs(r) >= 0 ) && ( abs(r) <= 1)
    rdelta = 2/3 - r^2 + 0.5 * abs(r)^3;
elseif (abs(r) > 1) && (abs(r) <= 2)
    rdelta = ( 2 - abs(r))^3 / 6;
else
    rdelta = 0;
end


end