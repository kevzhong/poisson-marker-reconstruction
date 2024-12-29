function [val,idx] = triangleTpt(data)
% Two endpoints on the curve "data"
x = [1 size(data,1)];
y = [data(1) data(end)];
% The slope of the line connecting the two endpoints
m = ( y(2) - y(1) )/( x(2) - x(1) );
pm= - 1 / m;
% Point on the curve (xc,yc), point on the line (xl,yl)
perpDist = zeros(size(data,1),1);
for i = 1:size(data,1)
    xc = i ; yc = data(i);
    yl = ( (m * xc) + (m^2 * yc) - (m * x(1)) + y(1) )/(1+ m^2);
    xl = xc - m*(yl - yc);
    % distance^2
    d2 = (xl - xc)^2 + (yl - yc)^2;
    perpDist(i) = d2;
end
[val, idx] = max(perpDist);

end