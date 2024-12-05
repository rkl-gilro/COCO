tests = convertCharsToStrings({'testing_cubes_Lab.mat', 'testing_spheres_Lab.mat', ...
'testing_adj1_Lab.mat','testing_adj3_Lab.mat','testing_matching_Lab.mat'});
ymax = 0;
for i=1:length(tests)
    aux = load(tests{i});
    results{i} = aux.CCI.cci_neutral(:, 2);
    participants_cell{i} = string(cell2mat(aux.Part_names));
    ymax = max(ymax, max(results{i}(:)));
end

%% Define colors: 1-blue, 2-yellow, 3-green, 4-red
K = 5;
c(:, 5, :) = repmat([30,144,255]./255, K, 1);
c(:, 3, :) = repmat([46,139,87]./255, K, 1);
c(:, 2, :) = repmat([218,165,32]./255, K, 1);
c(:, 4, :) = repmat([205,92,92]./255, K, 1);
c(:, 1, :) = repmat([103,58,183]./255, K, 1);

colors(:, :, 1) = c(:, :, 1);
colors(:, :, 2) = c(:, :, 2);
colors(:, :, 3) = c(:, :, 3);

figure;h = daviolinplot_color(results,'colors',colors,...
    'box',3,'boxcolor','w','scatter',2,'jitter',0,'scattercolor','same',...
'smoothing', .1,'boxspacing',1,'scattersize',220,'scatteralpha',0.7);
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]);
ylim([0, ymax+0.3]);
plot([xl(1)-0.1, xl(2)+0.2], [1 1], 'k--', 'LineWidth',3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
experiments = [];
experiments_num = [];
data = [];
participants = [];
for i=1:length(tests)
    experiments = [experiments;...
        repmat(convertCharsToStrings(tests{i}),length(results{i}),1)];                      % 
    experiments_num = [experiments_num;...
        repmat(i,length(results{i}),1)];
    aux2 = results{i};
    % aux2 = aux2(randperm(length(aux2)));
    data = [data;aux2];

    participants = [participants; participants_cell{i}];
end


%% We cannot perform anova because not all the participants run all the 
%% experimental conditions
%% RUN lme models instead: main effects and interactions
%% CCI ~ experiments + (1|participants)
%% BR shows differences between experiments but not CCI!!!
tbl = table(data,experiments, participants);
lme = fitlme(tbl,'data ~ experiments + (1|participants)');
lme2 = fitlme(tbl,'data ~ experiments + (experiments|participants)');
% lme3 = fitlme(tbl,'data ~ experiments + (experiments-1|participants)');
% lme4 = fitlme(tbl,'data ~ experiments -1 + (experiments-1|participants)');
% lme5 = fitlme(tbl,'data ~ experiments + (1|participants) + (experiments|participants)')

figure;
F = fitted(lme2);
R = response(lme2);
figure();
plot(R,F,'rx')
% xlabel('Response')
% ylabel('Fitted')
% figure();
% plotResiduals(lme2,'fitted')


each_part = unique(tbl.participants);

shapes = {'s', 'o', 'v', '^', 'hexagram'};
tbl = table(data,experiments, experiments_num, participants);
figure;
for i = 1:length(each_part)

    ind = find( tbl.participants == each_part(i)); 
    
    exp = tbl.experiments(ind);
    exp_num = tbl.experiments_num(ind);
    dat = tbl.data(ind);
    
    for j=1:length(exp_num)
    scatter(i, dat(j), 300, ...
        squeeze(c(1, exp_num(j), :))', 'filled', shapes{exp_num(j)});hold on
    end

end
xticks(1:length(each_part))
xticklabels(each_part)
axis([0 length(each_part)+1 0 1.6])
