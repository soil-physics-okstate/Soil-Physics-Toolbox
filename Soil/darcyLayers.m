function [Jw, pInterface, Keff] = darcyLayers (Ks, L, p, z, b)
%DARCYLAYERS function that calculates the flux per unit area of saoil with
%two or more different layers under staurated conditions.
%
% Ks = Saturated hydraulic conductivity [m/h].
% p = Hydrostatic potential [m].
% z = Gravitatinoal potential [m].
% b = Ponding depth [m].
%
% Example: [Jw,  Keff] = darcyLayers ([25 5], [75 25], [0 100], [0 100], 10);

%% Step 1. Calculate the hydraulic head (H).
H1 = z(:,1) + p(:,1);
H2 = z(:,2) + b;

%% Step 2. Calculate the effective saturated hydraulic conductivity (Keff).
Keff = sum(L) ./ sum(L./Ks);

%% Step 3. Apply Darcy's law.
Jw = -Keff .* (H2-H1)./(z(:,2)-z(:,1));

%% Step 4. Apply Darcy's law to each layer.
pInterface = Jw * (L-z(:,1)) ./ (-Ks) + L - H1;

end