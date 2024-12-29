function [rgb] = flattenColour(rgba)
%Covert 4 element rgba vector to a flattened 3 element rgb vector
%Default background colour set to white

opacity = rgba(4);
rgb = zeros(1,3);

rgb(1) = opacity * rgba(1) + (1 - opacity) * 1;
rgb(2) = opacity * rgba(2) + (1 - opacity) * 1;
rgb(3) = opacity * rgba(3) + (1 - opacity) * 1;


% const flattenedColor = {
%   r: opacity * fgCol.r + (1 - opacity) * bgCol.r,
%   g: opacity * fgCol.g + (1 - opacity) * bgCol.g,
%   b: opacity * fgCol.b + (1 - opacity) * bgCol.b,
% };
end

