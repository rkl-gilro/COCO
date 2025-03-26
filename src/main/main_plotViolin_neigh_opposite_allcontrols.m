clear
% close all

outdoor = 0;
relative = 0;
cond = 4;

if outdoor
    tests = ...
        convertCharsToStrings({'khaki_outdoor_ab.mat', 'rose_outdoor_ab.mat', ...
        'purple_outdoor_ab.mat','teal_outdoor_ab.mat',...
        'control_khaki_outdoor_ab.mat','control_rose_outdoor_ab.mat',...
        'control_purple_outdoor_ab.mat','control_teal_outdoor_ab.mat'});
else
    tests = ...
        convertCharsToStrings({'khaki_indoor_ab.mat', 'rose_indoor_ab.mat', ...
        'purple_indoor_ab.mat','teal_indoor_ab.mat',...
        'control_khaki_indoor_ab.mat','control_rose_indoor_ab.mat',...
        'control_purple_indoor_ab.mat','control_teal_indoor_ab.mat'});
end

colors_ill_clockwise = {'Yellow', 'Red','Blue', 'Green'};
colors_ill_anticlockwise = {'Green', 'Yellow', 'Red','Blue'};
colors_ill_oppclockwise = {'Blue', 'Green','Yellow', 'Red'};
colors_ill_oppanticlockwise = {'Red', 'Blue', 'Green','Yellow'};

order_illuminants = {'Neutral', 'Blue', 'Yellow','Green', 'Red'};

local = {"Khaki", "Rose", "Purple", "Teal", "Control"};

cci_surround = cell(cond, 1);
cci_control  = cell(cond, 1);
for i=1:cond%length(tests)
    cci_surround{i} = load(tests{i});
    cci_control{i}  = load(tests{cond+i});
end


cci_all = [];
participants_all = [];
experiment_all = [];
illuminants_all = [];
illuminants_all_rgb = [];
colorshifts_all = [];
neighbor_all = [];
local_control = [];

illuminants_names = {'Blue', 'Yellow', 'Green', 'Red'};
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
    aux = cci_surround{i}.CCI.ccip_neutral(:, 2:end);
    obs = cell2mat(cci_surround{i}.Part_names);

    obs_ = sum(str2num(obs(:, end-1:end)), 2)+1;

    data1{i} = aux(:);
    data1_all{i} = aux;
    data1_mean{i} = mean(aux, 2);

    aux2 = cci_control{i}.CCI.ccip_neutral(:, 2:end);
    data1_con{i} = aux2(:);
    data1_con_all{i} = aux2;
    data1_con_mean{i} = mean(aux2, 2);

    data1_all_con{i} = [aux(:) aux2(:)];
    
    if relative == 1
        data1_relative{i} = -1 + (data1{i})./(data1_con{i});
        data1_con_relative{i} = -1 + (data1_con{i})./(data1_con{i});

        aux3 = -1 + (aux)./(aux2);
        aux4 = -1 + (aux2)./(aux2);
    elseif relative == 2
        data1_relative{i} = 100.*(-1 + (data1{i})./(data1_con{i}));
        data1_con_relative{i} = 100.*(-1 + (data1_con{i})./(data1_con{i}));

        aux3 = 100.*(-1 + (aux)./(aux2));
        aux4 = 100.*(-1 + (aux2)./(aux2));
    else
        data1_relative{i} = -(data1_con{i}) + (data1{i});
        data1_con_relative{i} = -(data1_con{i}) + (data1_con{i});

        aux3 = -aux2 + (aux);
        aux4 = -aux2 + (aux2);
    end



    aux
    aux3
    % 
    % aux3 = -((aux))./aux2 + 1;
    % aux4 = -((aux2))./aux2 + 1;


    data1_rel_mean{i} = median(aux3, 2); 
    data1_con_rel_mean{i} = median(aux4, 2);

    data1_rel_mean_all{i} = [data1_con_rel_mean{i} data1_rel_mean{i}];
    

    data1_rel_meanv(obs_,i) = median(aux3, 2); 
    data1_con_rel_meanv(obs_,i) = median(aux4, 2);

    data1_rel_meanv(data1_rel_meanv == 0) = NaN;
    data1_rel_meanv_cell{i} = data1_rel_meanv(:,i);

end

% for i=1:length(cci_surround)
%     data1_relative{i} = 1 - abs(data1{i} - mean(data1_con{i}))./ ...
%         mean(data1_con{i});
%     data1_con_relative{i} = 1 - abs(data1_con{i} - mean(data1_con{i}))./ ...
%         mean(data1_con{i});
% 
% 
%     data1_rel_mean{i} = mean(aux, 2);
%     data1_con_rel_mean{i} = mean(aux2, 2);
% end
figure;
plot(0:6, zeros(1, 7), 'k--','LineWidth',3);hold on
h = daboxplot(data1_rel_meanv_cell,'xtlabels', condition_names,'whiskers',0,...
'scatter',1,'scattersize',250,'scatteralpha',0.6,'outliers',0,...
'colors', c,'linkline',1, 'withinline', 1); % experiments_ill.^.45
set(gca, 'FontSize', 40, 'fontname','L M Roman10');
ylim([-0.4 .4])
yticks([-.5 0 1.2])
ylabel('Absolute error (\Delta AE)')

figure;
plot(0:6, zeros(1, 7), 'k--','LineWidth',3);hold on
h = dabarplot(data1_relative,'xtlabels', condition_names,...
'colors',c,'errorbars','IQR','numbers',1,'round',3, 'scattersize', 40);
set(gca, 'FontSize', 40, 'fontname','L M Roman10');
ylim([-0.3 .2])
yticks([-.3 -.2 -.1 0 .1 .2])

aa = data1_con_rel_meanv;
aa(aa==0) = NaN;
bb = data1_rel_meanv;
bb(bb==0) = NaN;

force = [median(aa,2,'omitnan') ...
    median(bb,2,'omitnan')];
figure;
h = daboxplot(force,'xtlabels', {'Control','LocalSurround'},'whiskers',0,...
'scatter',1,'scattersize',250,'scatteralpha',0.6,'outliers',0,...
'colors',[.7 .7 .6],'linkline',1,'withinlines', 1);
set(gca, 'FontSize', 30, 'fontname','L M Roman10');
ylim([-0.6 .5])

figure;h = daviolinplot(data1_rel_mean,'colors',c,'box',3,...
'boxcolor','w','scatter',2,'jitter',0,'scattercolor','same',...
'smoothing', .1,'boxspacing',1.3, ...
'scattersize',500,'scatteralpha',0.7,...
'xtlabels', condition_names,'linkline',1,'withinlines', 1);
ylim([0.5 1.25])
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]);
set(gca, 'FontSize', 50, 'fontname','L M Roman10');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define table for statistical analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(cci_surround)
    aux = [cci_surround{i}.CCI.ccip_neutral(:, 2:end); ...
           cci_control{i}.CCI.ccip_neutral(:, 2:end)];


    aux2 = [(cci_surround{i}.Part_names);
            cci_control{i}.Part_names];
    
    size_aux = cci_surround{i}.CCI.ccip_neutral(:, 2:end);
    aux_part = size(size_aux, 1);

    size_aux2 = cci_control{i}.CCI.ccip_neutral(:, 2:end);
    aux_part2 = size(size_aux2, 1);


    aux_ill = [repmat("Blue",aux_part,1);repmat("Yellow",aux_part,1);...
        repmat("Green",aux_part,1);repmat("Red",aux_part,1);...
        repmat("Blue",aux_part2,1);repmat("Yellow",aux_part2,1);...
        repmat("Green",aux_part2,1);repmat("Red",aux_part2,1)];
    
    local_control = [local_control;
        repmat("Local Surround",length(size_aux(:)),1);
        repmat("Control",length(size_aux2(:)),1)];

    cci_all = [cci_all; size_aux(:); size_aux2(:)];
    illuminants_all = [illuminants_all; aux_ill];

    experiment_all = [experiment_all; 
                      repmat(local{i},length(size_aux(:)),1);
                      repmat(strcat('Control',local{i}),...
                      length(size_aux2(:)),1)];

    aux_partcell = repmat(aux2, 1,size(aux,2));
    participants_all = [participants_all; aux_partcell(:)];

    switch i
        case 1 % Khaki
            aux_neig = [repmat("Opposite",aux_part,1);...
                repmat("Neighbor",aux_part,1);...
        repmat("Neighbor",aux_part,1);repmat("Opposite",aux_part,1);
        repmat("NoneO",aux_part2,1);...
                repmat("NoneN",aux_part2,1);...
        repmat("NoneN",aux_part2,1);repmat("NoneO",aux_part2,1);];

        case 2 % Rose
            aux_neig = [repmat("Opposite",aux_part,1);...
                repmat("Neighbor",aux_part,1);...
        repmat("Opposite",aux_part,1);repmat("Neighbor",aux_part,1);
        repmat("NoneO",aux_part2,1);...
                repmat("NoneN",aux_part2,1);...
        repmat("NoneO",aux_part2,1);repmat("NoneN",aux_part2,1);];

        case 3 % Purple
            aux_neig = [repmat("Neighbor",aux_part,1);...
                repmat("Opposite",aux_part,1);...
        repmat("Opposite",aux_part,1);repmat("Neighbor",aux_part,1);
        repmat("NoneN",aux_part2,1);...
                repmat("NoneO",aux_part2,1);...
        repmat("NoneO",aux_part2,1);repmat("NoneN",aux_part2,1);];

        case 4 % Teal
            aux_neig = [repmat("Neighbor",aux_part,1);...
                repmat("Opposite",aux_part,1);...
        repmat("Neighbor",aux_part,1);repmat("Opposite",aux_part,1);
        repmat("NoneN",aux_part2,1);...
                repmat("NoneO",aux_part2,1);...
        repmat("NoneN",aux_part2,1);repmat("NoneO",aux_part2,1);];

    end
    neighbor_all = [neighbor_all; aux_neig(:)];
end

tablenames = ...
    {'CCI', 'Participant', 'LocalSurround','Baseline', 'Illuminant', 'Direction'};
tbl = table(cci_all, ...
    participants_all, experiment_all, local_control, ...
    illuminants_all, neighbor_all,'VariableNames', tablenames);

tbl.LocalSurround = nominal(tbl.LocalSurround);
tbl.Illuminant    = nominal(tbl.Illuminant);
tbl.Direction     = nominal(tbl.Direction);
tbl.Baseline      = nominal(tbl.Baseline);

%% Statistical analysis based on the direction: Neighboring, Opposing
lme{1} = fitlme(tbl, 'CCI ~ 1 + Direction + (1|Participant)');  
lme{2} = fitlme(tbl, 'CCI ~ 1 + Baseline + (1|Participant)');  

glme{1} = fitglme(tbl, 'CCI ~ 1 + Direction + (1|Participant)');
glme{2} = fitglme(tbl, 'CCI ~ 1 + Baseline + (1|Participant)'); 

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
    '(1|Participant)']);  
lme_ls{2} = fitlme(tbl, ['CCI ~ 1 + LocalSurround*Illuminant + ' ...
    '(1|Participant)']); 
lme_ls{3} = fitlme(tbl, ['CCI ~ 1 + Baseline*Illuminant + ' ...
    'Baseline^2 + (1|Participant)']); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot CCI results based on the neighboring, opposing illuminants
c2 = [125,125,155;
     176,196,222;
     112,128,144;
     125,125,125;
     176,196,222]./255;

data2{1} = tbl.CCI(tbl.Direction == 'NoneO');
data2{2} = tbl.CCI(tbl.Direction == 'NoneN');
data2{3} = tbl.CCI(tbl.Direction == 'Opposite');
data2{4} = tbl.CCI(tbl.Direction == 'Neighbor');

data_none = tbl.CCI(tbl.Direction == 'NoneO' | tbl.Direction == 'NoneN');

if relative
    data3{1} = -1 + (data2{1}) ./ data2{1};
    data3{2} = -1 + (data2{2}) ./ data2{2};
    data3{3} = -1 + (data2{3}) ./ data2{1}; % ./ mean(data2{1})
    data3{4} = -1 + (data2{4}) ./ data2{2};
else
    data3{1} = data2{1} - data2{1};
    data3{2} = data2{2} - data2{2};
    data3{3} = data2{3} - data2{1}; % ./ mean(data2{1})
    data3{4} = data2{4} - data2{2};
end

unique_part = convertCharsToStrings(unique(participants_all));

for i=1:length(unique_part)
    
    opp(i) = mean(tbl.CCI(tbl.Direction == 'Opposite' & ...
        convertCharsToStrings(tbl.Participant) == unique_part(i, :)));
    neig(i) = mean(tbl.CCI(tbl.Direction == 'Neighbor' & ...
        convertCharsToStrings(tbl.Participant) == unique_part(i, :)));

    dnoneo(i) = mean(tbl.CCI(tbl.Direction == 'NoneO' & ...
        convertCharsToStrings(tbl.Participant) == unique_part(i, :)));

    dnonen(i) = mean(tbl.CCI(tbl.Direction == 'NoneN' & ...
        convertCharsToStrings(tbl.Participant) == unique_part(i, :)));

end

all_ls = nan(length(unique_part), length(local)-1, ...
    length(illuminants_names));
all_c_ls = nan(length(unique_part), length(local)-1, ...
    length(illuminants_names));
for i=1:length(unique_part)

    for j=1:length(local)-1
        ls = convertCharsToStrings(local{j});
        c_ls = strcat('Control', ls);

        for k=1:length(illuminants_names)
            if ~isempty(tbl.CCI(tbl.LocalSurround == ls & ...
                    convertCharsToStrings(tbl.Participant) == unique_part(i, :) & ...
                    tbl.Illuminant == illuminants_names(k)))

                all_ls(i, j, k) = tbl.CCI(tbl.LocalSurround == ls & ...
                    convertCharsToStrings(tbl.Participant) == unique_part(i, :) & ...
                    tbl.Illuminant == illuminants_names(k));
            end

            if ~isempty(tbl.CCI(tbl.LocalSurround == c_ls & ...
                    convertCharsToStrings(tbl.Participant) == unique_part(i, :) & ...
                    tbl.Illuminant == illuminants_names(k)))

                all_c_ls(i, j, k) = tbl.CCI(tbl.LocalSurround == c_ls & ...
                    convertCharsToStrings(tbl.Participant) == unique_part(i, :) & ...
                    tbl.Illuminant == illuminants_names(k));
            end

        end
        
    end

end


figure;
plot(0:6, ones(1, 7), 'k--','LineWidth',3);hold on
h = daboxplot(squeeze(all_ls(:, 1, :)),...
    'xtlabels', illuminants_names,'whiskers',0,...
'scatter',1,'scattersize',350,'scatteralpha',0.6,'withinlines',1,...
'linkline',1,'outliers',0,'colors',c);
set(gca, 'FontSize', 30, 'fontname','L M Roman10');

figure;
count = 1;
for i=1:length(illuminants_names)
    subplot(2, 2, count)
    plot(0:6, ones(1, 7), 'k--','LineWidth',3);hold on
    h = daboxplot(squeeze(all_ls(:, :, i)),...
        'xtlabels', local,'whiskers',0,...
        'scatter',1,'scattersize',350,'scatteralpha',0.6,'withinlines',1,...
        'linkline',1,'outliers',0,...
        'colors',(2.*experiments_ill(i, :)).^.45);
    ylim([0.0 1.1])
    set(gca, 'FontSize', 30, 'fontname','L M Roman10');
    count = count + 1;
end

figure;
count = 1;
for i=1:length(illuminants_names)
    subplot(2, 2, count)
    plot(0:6, ones(1, 7), 'k--','LineWidth',3);hold on
    h = daboxplot(squeeze(all_c_ls(:, :, i)),...
        'xtlabels', local,'whiskers',0,...
        'scatter',1,'scattersize',350,'scatteralpha',0.6,'withinlines',1,...
        'linkline',1,'outliers',0,...
        'colors',(2.*experiments_ill(i, :)).^.45);
    ylim([0.0 1.1])
    set(gca, 'FontSize', 30, 'fontname','L M Roman10');
    count = count + 1;
end

data5(:,1) = data3{3};
data5(:,2) = data3{4};

if relative
    data7(:,1) = -1 + (opp)./dnoneo;
    data7(:,2) = -1 + (neig)./dnonen;
else
    data7(:,1) = -dnoneo + (opp);
    data7(:,2) = -dnonen + (neig);
end

figure;
plot(0:6, zeros(1, 7), 'k--','LineWidth',3);hold on
h = daboxplot(data7,'xtlabels', {"Opposing","Neighboring"},'whiskers',0,...
'scatter',1,'scattersize',350,'scatteralpha',0.6,'withinlines',1,...
'linkline',1,'outliers',0,'colors',c2);
set(gca, 'FontSize', 30, 'fontname','L M Roman10');
ylim([-.6 .6])
ylabel('Absolute error (\Delta AE)')

figure;h2 = daviolinplot(data2,'colors',c2,'box',1,...
'boxcolor','w','scatter',2,'jitter',1,'scattercolor','same',...
'smoothing', .15,'boxspacing',1.3, ...
'scattersize',200,'scatteralpha',0.7,'linkline',1,...
'xtlabels', {"ControlO","ControlN","Opposing","Neighboring"}); % condition_names2
ylim([0 1.2])

figure;h2 = daviolinplot(data5,'colors',c,'box',1,...
'boxcolor','w','scatter',2,'jitter',1,'scattercolor','same',...
'smoothing', .11,'boxspacing',1.3, ...
'scattersize',200,'scatteralpha',0.7,'linkline',1,...
'xtlabels', {"Opposing","Neighboring"}); % condition_names2
ylim([0 1.2])
set(gca, 'FontSize', 20, 'fontname','L M Roman10');

figure;h2 = daviolinplot(data3,'colors',c2,'box',1,...
'boxcolor','w','scatter',2,'jitter',1,'scattercolor','same',...
'smoothing', .15,'boxspacing',1.3, ...
'scattersize',200,'scatteralpha',0.7,'linkline',1,...
'xtlabels', {"ControlO","ControlN","Opposing","Neighboring"}); % condition_names2
ylim([0 1.5])
% ylim([-150 100])
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]);
set(gca, 'FontSize', 20, 'fontname','L M Roman10');

figure;
h = daboxplot(data3,'xtlabels', ...
    {"ControlO","ControlN","Opposing","Neighboring"},'whiskers',0,...
    'scatter',1,'scattersize',250,'scatteralpha',0.6,'withinlines',...
    1,'outliers',0,'colors',c2);

figure;
h = daviolinplot(data5,'colors',c2,'boxcolors','k','outliers',0,...
    'box',3,'boxwidth',0.8,'scatter',2,'scattersize',35,'jitter',1,...
    'xtlabels', {"ControlO","ControlN","Opposing","Neighboring"}); 