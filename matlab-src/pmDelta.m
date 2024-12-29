function rdelta = pmDelta(r)
% Compute the Peksin & McQueen (1989) delta function

d1 = @(r) 3 - 2 * abs(r) + sqrt(   1 + 4 * abs(r) - 4 * r^2   );

if abs(r) < 1
    rdelta = d1(r) / 8;
elseif (abs(r) > 1) && (abs(r) < 2)
    rdelta = 0.5 - d1(2 - abs(r));
else
    rdelta = 0;
end



end