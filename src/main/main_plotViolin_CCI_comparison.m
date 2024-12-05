%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Plot comparison between different CCIs
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tests = {'testing_matching_Lab.mat'};
ymax = 0;
labels = {'CCI','CCI neutral','BR','BR neutral'};

aux = load(tests{1});
cci_table = struct2table(aux.CCI);
for i=1:size(cci_table,2)

    aux2 = table2array(cci_table(:, i));
    results{i} = aux2(:, 2);

    ymax = max(ymax, max(results{i}(:)));
end

%% Define colors: 1-blue, 2-yellow, 3-green, 4-red
K = 1;
c(1, :) = repmat([229, 115, 115]./255, K, 1);
c(2, :) = repmat([239, 83, 80]./255, K, 1);
c(3, :) = repmat([77, 208, 225]./255, K, 1);
c(4, :) = repmat([0, 151, 167]./255, K, 1);

figure;h = daviolinplot(results,'colors',c,'box',3,...
'boxcolor','w','scatter',2,'jitter',0,'scattercolor','same',...
'smoothing', .15,'boxspacing',1, ...
'scattersize',320,'scatteralpha',0.7);
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]);
ylim([0, 1.5]);
plot([xl(1)-0.1, xl(2)+0.2], [1 1], 'k--', 'LineWidth',3);
xticklabels(labels)
set(gca, 'FontSize', 20, 'fontname','Times New Roman');