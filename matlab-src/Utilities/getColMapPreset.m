function Cmap = getColMapPreset(Pr,nlevels)
%Custom colormap presets

%RGB = white to:

%[1,0.0551181102362205,0] Pr050
%[0,0,0] %Something else
%[0,0.0393700787401575,0.968253968253968] %Pr100
%[0,0.425196850393701,0] %Pr200

% Source => Target = (BGColor + Source) =
% Target.R = ((1 - Source.A) * BGColor.R) + (Source.A * Source.R)
% Target.G = ((1 - Source.A) * BGColor.G) + (Source.A * Source.G)
% Target.B = ((1 - Source.A) * BGColor.B) + (Source.A * Source.B)

startcol = [ 1 1 1];

nlevels = nlevels + 1;
switch Pr
    case 0.5
        endcol = [1,0.0551181102362205,0];
    case 1.0
        endcol = [0,0.0393700787401575,0.968253968253968];
    case 2.0
        endcol = [0,0.425196850393701,0];
    otherwise
        endcol = [0 0 0];
end


Cmap = [linspace(startcol(1),endcol(1),nlevels)', linspace(startcol(2),endcol(2),nlevels)',...
    linspace(startcol(3),endcol(3),nlevels)'];

Cmap = Cmap(2:end,:);

% opacities = linspace(0.2,1,nlevels);
% Cmap = [repmat(endcol,[nlevels,1]), opacities'];

end

