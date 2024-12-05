tests = {'control_khaki_indoor.mat','khaki_indoor.mat'};
ymax = 0;
for i=1:length(tests)
    aux = load(tests{i});
    results{i} = aux.CCI.cci_neutral(:, 2:end);

    ymax = max(ymax, max(results{i}(:)));
end

%% Define colors: 1-blue, 2-yellow, 3-green, 4-red
K = 4;
c(:, 1, :) = repmat([30,144,255]./255, K, 1);
c(:, 3, :) = repmat([46,139,87]./255, K, 1);
c(:, 2, :) = repmat([218,165,32]./255, K, 1);
c(:, 4, :) = repmat([205,92,92]./255, K, 1);

colors(:, :, 1) = c(:, :, 1);
colors(:, :, 2) = c(:, :, 2);
colors(:, :, 3) = c(:, :, 3);

figure;h = daviolinplot_color(results,'colors',colors,'box',3,...
'boxcolor','w','scatter',2,'jitter',0,'scattercolor','same',...
'smoothing', .1,'boxspacing',1, ...
'scattersize',220,'scatteralpha',0.7);
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]);
ylim([0, 1.2]);
plot([xl(1)-0.1, xl(2)+0.2], [1 1], 'k--', 'LineWidth',3);