function [ddzbot,ddztop] = get_ddz_walls(dz_ghost, Ug)

n2 = size(Ug,1);
if n2 - length(dz_ghost) == 6
    dz_ghost = get_dz_ghost(dz_ghost);
end

a1 = 9/8;
a2 = -1/24;

ddzbot = ( a1 * (Ug(4)-Ug(3)) + ...
             a2 * (Ug(5)-Ug(2)) ) / ...
           (0.5*dz_ghost(4) + 0.5*dz_ghost(3));

%if isfull_dz(dz_ghost)
    Ug = flipud(Ug);
    ddztop = ( a1 * (Ug(4)-Ug(3)) + ...
                 a2 * (Ug(5)-Ug(2)) ) / ...
               (0.5*dz_ghost(4) + 0.5*dz_ghost(3));
%end

%ddz = mean(ddz);
end
