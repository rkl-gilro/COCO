function [x0, y0] = daylightlocus_T2xy_rec(T0, show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Main compute xy chromaticities from daylight locus 
%% given its T temperature
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% xy Daylight locus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    show = 0;
end
%% Rectify T0
T0 = T0 *1.438776877/1.438;

T = 4000:25000;
xc = zeros(1, length(T) );
yc = zeros(1, length(T) );

Taux = T(T <= 7000);
xc(1:find(T==7000)) = -4.607.*(10^9./Taux.^3) + 2.9678.*(10^6./Taux.^2) ...
                     + 0.09911.*(10^3./Taux) + 0.244063;
clear Taux
Taux = T(T > 7000);
xc(find(T==7000)+1:end) = -2.0064.*(10^9./Taux.^3) + 1.9018.*(10^6./Taux.^2) ...
                     + 0.24748.*(10^3./Taux) + 0.237040;

yc = -3.*(xc.^2) + 2.87.*xc - 0.275;

x0 = mean(xc(T >= floor(T0) & T <=round(T0))); 
y0 = mean(yc(T >= floor(T0) & T <=round(T0)));


if show
    figure;
    plotChromaticity();hold on;
    plot(xc, yc, 'k-', 'LineWidth',2);hold on;
    plot(xc(T == 6500), yc(T == 6500), 'b+', 'LineWidth',2);
    plot(x0, y0, 'ko');

    set(gca,  'FontSize', 20, 'fontname','Times New Roman');
    grid on
    set(gcf,'renderer','Painters');
end