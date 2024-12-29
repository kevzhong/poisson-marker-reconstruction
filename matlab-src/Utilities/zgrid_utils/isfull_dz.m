function bool = isfull_dz(dz)
% check if dz is more or less symmetric
bool = 0;

if abs(find(dz==max(dz),1) - length(dz)/2) < 2 % central maximum
    bool = 1;
elseif dz(round(end/2)+1) == dz(round(end/2)) % mirrored
    bool = 1;
end

end