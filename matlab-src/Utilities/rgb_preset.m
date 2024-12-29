function [rgb] = rgb_preset()
%Get preset colors (for manuscript)

fun = @(m)Lab_to_DIN99(srgb_to_Lab(m)); % DIN99 = good colorspace (close colors).
[rgb,~] = maxdistcolor(4,fun);

%Switch 2 and 4

rgb([2 4],:) = rgb([4 2],:);
end

