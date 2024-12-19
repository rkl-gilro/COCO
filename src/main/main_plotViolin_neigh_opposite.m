clear
% close all
tests = ...
    convertCharsToStrings({'khaki_indoor_ab.mat', 'rose_indoor_ab.mat', ...
'purple_indoor_ab.mat','teal_indoor_ab.mat',...
'control_khaki_indoor_ab.mat'});
% tests = ...
%     convertCharsToStrings({'khaki_outdoor_ab.mat', 'rose_outdoor_ab.mat', ...
% 'purple_outdoor_ab.mat','teal_outdoor_ab.mat',...
% 'control_khaki_outdoor_ab.mat'});

outdoor = 1;
relative = 1;

colors_ill_clockwise = {'Yellow', 'Red','Blue', 'Green'};
colors_ill_anticlockwise = {'Green', 'Yellow', 'Red','Blue'};
colors_ill_oppclockwise = {'Blue', 'Green','Yellow', 'Red'};
colors_ill_oppanticlockwise = {'Red', 'Blue', 'Green','Yellow'};

order_illuminants = {'Neutral', 'Blue', 'Yellow','Green', 'Red'};

local = {"Khaki", "Rose", "Purple", "Teal", "Control"};

cci_surround = cell(length(tests), 1);
for i=1:length(tests)
    cci_surround{i} = load(tests{i});
end


cci_all = [];
participants_all = [];
experiment_all = [];
illuminants_all = [];
illuminants_all_rgb = [];
colorshifts_all = [];
neighbor_all = [];

experiments_ill = [0.2316 0.2833 0.4579; %blue
    0.38 0.2441 0.1305;  %yellow
    0.233 0.3 0.2003;  %green
    0.4152 0.2114 0.3961]; %red
%    0.2912 0.2705 0.2614; %neutral

%% Plot the local surround CCI and the control %%%%%%%%%%%%%%%%%%%%%%%%%%%%
condition_names = local;
condition_names2 = {"Control", "Opposing", "Neighboring"};

c = [189,183,107;
    205,92,92;
    147,112,219;
    95,158,160;
    125,125,125;]./255;

for i=1:length(cci_surround)
    % data1{i} = mean(cci_surround{i}.CCI.cci_neutral(:, 2:end)', 2);
    aux = cci_surround{i}.CCI.cci_neutral(:, 2:end);
    data1{i} = aux(:);
end

figure;h = daviolinplot(data1,'colors',c,'box',3,...
'boxcolor','w','scatter',2,'jitter',0,'scattercolor','same',...
'smoothing', .1,'boxspacing',1.3, ...
'scattersize',200,'scatteralpha',0.7,...
'xtlabels', condition_names);
ylim([0 1.2])
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]);
set(gca, 'FontSize', 20, 'fontname','Times New Roman');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define table for statistical analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(cci_surround)
    aux = cci_surround{i}.CCI.cci_neutral(:, 2:end);
    aux2 = cci_surround{i}.Part_names;
    aux_part = size(aux, 1);
    aux_ill = [repmat("Blue",aux_part,1);repmat("Yellow",aux_part,1);...
        repmat("Green",aux_part,1);repmat("Red",aux_part,1)];

    cci_all = [cci_all; aux(:)];
    illuminants_all = [illuminants_all; aux_ill];

    experiment_all = [experiment_all; repmat(local{i},length(aux(:)),1)];

    aux_partcell = repmat(aux2, 1,size(aux,2));
    participants_all = [participants_all; aux_partcell(:)];

    switch i
        case 1 % Khaki
            aux_neig = [repmat("Opposite",aux_part,1);...
                repmat("Neighbor",aux_part,1);...
        repmat("Neighbor",aux_part,1);repmat("Opposite",aux_part,1)];

        case 2 % Rose
            aux_neig = [repmat("Opposite",aux_part,1);...
                repmat("Neighbor",aux_part,1);...
        repmat("Opposite",aux_part,1);repmat("Neighbor",aux_part,1)];

        case 3 % Purple
            aux_neig = [repmat("Neighbor",aux_part,1);...
                repmat("Opposite",aux_part,1);...
        repmat("Opposite",aux_part,1);repmat("Neighbor",aux_part,1)];

        case 4 % Teal
            aux_neig = [repmat("Neighbor",aux_part,1);...
                repmat("Opposite",aux_part,1);...
        repmat("Neighbor",aux_part,1);repmat("Opposite",aux_part,1)];

        case 5
            aux_neig = [repmat("None",aux_part,1);...
                repmat("None",aux_part,1);...
        repmat("None",aux_part,1);repmat("None",aux_part,1)];
    end
    neighbor_all = [neighbor_all; aux_neig];
end

tablenames = ...
    {'CCI', 'Participant', 'LocalSurround', 'Illuminant', 'Direction'};
tbl = table(cci_all, ...
    participants_all, experiment_all, ...
    illuminants_all, neighbor_all,'VariableNames', tablenames);

tbl.LocalSurround = nominal(tbl.LocalSurround);
tbl.Illuminant    = nominal(tbl.Illuminant);
tbl.Direction     = nominal(tbl.Direction);

%% Statistical analysis based on the direction: Neighboring, Opposing
lme{1} = fitlme(tbl, 'CCI ~ 1 + Direction + (Direction|Participant)');
lme{2} = fitlme(tbl, 'CCI ~ 1 + Direction + (1|Participant)');
lme{3} = fitlme(tbl, 'CCI ~ 1 + Direction + (Direction-1|Participant)');
lme{4} = fitlme(tbl, ['CCI ~ 1 + Direction + (1 + LocalSurround|Participant)' ...
    '+ (1 + Illuminant|Participant)']);
lme{5} = fitlme(tbl, ['CCI ~ 1 + Direction + (1 + LocalSurround|Participant)' ...
    '+ (1 + Illuminant|Participant) + (1+Direction|Participant)']);

glme{1} = fitglme(tbl, 'CCI ~ 1 + Direction + (Direction|Participant)');
glme{2} = fitglme(tbl, 'CCI ~ 1 + Direction + (1|Participant)');
glme{3} = fitglme(tbl, 'CCI ~ 1 + Direction + (Direction-1|Participant)');
glme{4} = fitglme(tbl, ['CCI ~ 1 + Direction + (1 + LocalSurround|Participant)' ...
'+ (1 + Illuminant|Participant)'], 'FitMethod', 'Laplace');
glme{5} = fitglme(tbl, ['CCI ~ 1 + Direction + (1 + LocalSurround|Participant)' ...
'+ (1 + Illuminant|Participant) + (1+Direction|Participant)'], ...
'FitMethod', 'Laplace');

%% Plot the fitting of the models in a figure
shapes = {'o', 'v'};

sz = 100;
figure;
plot([0 1], [0 1], 'k--');hold on
for i=1:length(lme)
    F = fitted(lme{i});
    R = response(lme{i});
    scatter(R,F,sz,'MarkerFaceColor', c(i, :), 'Marker', shapes{1});hold on
    
    corrcoef(F,R)
    

    F = fitted(glme{i});
    R = response(glme{i});
    scatter(R,F,sz,'Marker',shapes{2},'MarkerFaceColor', c(i, :));hold on
    
    corrcoef(F,R)
end
axis([0 1 0 1])
axis square

% Uncomment for comparing models
% compare(lme, lme2)

%% Statistical analysis based on the local surround colour: 
%% 
lme_ls{1} = fitlme(tbl, ['CCI ~ 1 + LocalSurround + ' ...
    '(LocalSurround|Participant)']);
lme_ls{2} = fitlme(tbl, ['CCI ~ 1 + LocalSurround + ' ...
    '(1|Participant)']);
lme_ls{3} = fitlme(tbl, ['CCI ~ 1 + LocalSurround + ' ...
    '(LocalSurround-1|Participant)']);
lme_ls{4} = fitlme(tbl, ['CCI ~ 1 + LocalSurround + ' ...
    '(1 + LocalSurround|Participant)' ...
    '+ (1 + Illuminant|Participant)']);
lme_ls{5} = fitlme(tbl, ['CCI ~ 1 + LocalSurround + ' ...
    '(1 + LocalSurround|Participant)' ...
    '+ (1 + Illuminant|Participant) + (1+Direction|Participant)']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot CCI results based on the neighboring, opposing illuminants
c2 = [125,125,125;
     176,196,222;
     112,128,144]./255;

data2{1} = tbl.CCI(tbl.Direction == 'None');
data2{2} = tbl.CCI(tbl.Direction == 'Opposite');
data2{3} = tbl.CCI(tbl.Direction == 'Neighbor');
figure;h2 = daviolinplot(data2,'colors',c2,'box',3,...
'boxcolor','w','scatter',2,'jitter',0,'scattercolor','same',...
'smoothing', .1,'boxspacing',1.3, ...
'scattersize',200,'scatteralpha',0.7,...
'xtlabels', condition_names2);
ylim([0 1.2])
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]);
set(gca, 'FontSize', 20, 'fontname','Times New Roman');