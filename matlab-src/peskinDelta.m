function rdelta = peskinDelta(r)
% Compute the Peksin (1977) delta function


if abs(r) < 2
    rdelta = 0.25 * ( 1 + cos(pi*r/2) );
else
    rdelta = 0;
end