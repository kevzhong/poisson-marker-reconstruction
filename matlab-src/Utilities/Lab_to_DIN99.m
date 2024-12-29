function Lab99 = Lab_to_DIN99(Lab)
% Convert a matrix of L*a*b* values to DIN99 values (DIN 6176).
%
% (c) 2018-2019 Stephen Cobeldick
%
%%% Syntax:
% DIN99 = Lab_to_DIN99(Lab)
%
%% Inputs and Outputs
%
%%% Input Argument:
% Lab = Numeric Matrix, size Nx3, with CIE L*a*b* values [L,a,b,].
%
%%% Output Argument:
% Lab99 = Numeric Matrix, size Nx3, with DIN99 values [L99,a99,b99].
%
% See also SRGB_TO_LAB SRGB_TO_JAB MAXDISTCOLOR MAXDISTCOLOR_VIEW

L99 = 105.51 * log(1 + 0.0158*Lab(:,1));
e =     (Lab(:,2).*cosd(16)+Lab(:,3).*sind(16));
f = 0.7*(Lab(:,3).*cosd(16)+Lab(:,2).*sind(16));
G = sqrt(e.^2 + f.^2);
C99 = log(1 + 0.045*G)./0.045;
h99 = atan2(f,e);
a99 = C99 .* cos(h99);
b99 = C99 .* sin(h99);
Lab99 = [L99,a99,b99];
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Lab_to_DIN99