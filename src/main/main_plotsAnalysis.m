clear
% close all

outdoor = 0;
relative = 0;
cond = 4;

cci_method = 4; % 1-ccip, 2-ccip neutral, 3-BRp, 4-BRp neutral
cci_method_name = {'ccip', 'ccipn', 'BRp', 'BRpn', 'karl',...
                   'cci', 'ccin', 'BR', 'BRn'};

if outdoor
    save_table_filename = strcat('LocalSurround_Control_outdoor_', ...
    cci_method_name{cci_method});
else
    save_table_filename = strcat('LocalSurround_Control_indoor_', ...
    cci_method_name{cci_method});
end

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

colors_rgb_ls = [0.2331    0.2311    0.1376;
    0.4340    0.1623    0.1963;
    0.2630    0.1972    0.3911;
    0.1264    0.2485    0.3046];

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
colorshifts_control_all = [];
colorshifts_localsurround_all = [];
neighbor_all = [];
space_competitors = [];
local_control = [];
scene = [];
diff = [];

illuminants_names = {'Blue', 'Yellow', 'Green', 'Red'};
experiments_ill = [0.2912 0.2705 0.2614; %neutral
    0.2316 0.2833 0.4579; %blue
    0.38 0.2441 0.1305;  %yellow
    0.233 0.3 0.2003;  %green
    0.4152 0.2114 0.3961]; %red
%    ; %neutral

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

    switch cci_method
        case 1
            aux = cci_surround{i}.CCI.ccip(:, 2:end);
        case 2
            aux = cci_surround{i}.CCI.ccip_neutral(:, 2:end);
        case 3
            aux = cci_surround{i}.CCI.BRp(:, 2:end);
        case 4
            aux = cci_surround{i}.CCI.BRp_neutral(:, 2:end);
        case 5
            aux = cci_surround{i}.CCI.Karl(:, 2:end);
        case 6
            aux = cci_surround{i}.CCI.cci(:, 2:end);
        case 7
            aux = cci_surround{i}.CCI.cci_neutral(:, 2:end);
        case 8
            aux = cci_surround{i}.CCI.BR(:, 2:end);
        case 9
            aux = cci_surround{i}.CCI.BR_neutral(:, 2:end);
    end

    % aux = cci_surround{i}.CCI.ccip_neutral(:, 2:end);
    obs = cell2mat(cci_surround{i}.Part_names);

    obs_ = sum(str2num(obs(:, end-1:end)), 2)+1;

    data1{i} = aux(:);
    data1_all{i} = aux;
    data1_mean{i} = mean(aux, 2);
    
    switch cci_method
        case 1
            aux2 = cci_control{i}.CCI.ccip(:, 2:end);
        case 2
            aux2 = cci_control{i}.CCI.ccip_neutral(:, 2:end);
        case 3
            aux2 = cci_control{i}.CCI.BRp(:, 2:end);
        case 4
            aux2 = cci_control{i}.CCI.BRp_neutral(:, 2:end);
        case 5
            aux2 = cci_control{i}.CCI.Karl(:, 2:end);
        case 6
            aux2 = cci_control{i}.CCI.cci(:, 2:end);
        case 7
            aux2 = cci_control{i}.CCI.cci_neutral(:, 2:end);
        case 8
            aux2 = cci_control{i}.CCI.BR(:, 2:end);
        case 9
            aux2 = cci_control{i}.CCI.BR_neutral(:, 2:end);
    end

    % aux2 = cci_control{i}.CCI.ccip_neutral(:, 2:end);
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


    data1_rel_mean{i} = median(aux3, 2); 
    data1_con_rel_mean{i} = median(aux4, 2);

    data1_rel_mean_all{i} = [data1_con_rel_mean{i} data1_rel_mean{i}];
    

    data1_rel_meanv(obs_,i) = median(aux3, 2); 
    data1_con_rel_meanv(obs_,i) = median(aux4, 2);

    data1_rel_meanv(data1_rel_meanv == 0) = NaN;
    data1_rel_meanv_cell{i} = data1_rel_meanv(:,i);

end


figure;
plot(0:6, zeros(1, 7), 'k--','LineWidth',3);hold on
h = daboxplot(data1_rel_meanv_cell,'xtlabels', condition_names,'whiskers',0,...
'scatter',1,'scattersize',250,'scatteralpha',0.6,'outliers',0,...
'colors', c,'linkline',1, 'withinline', 1); % experiments_ill.^.45
set(gca, 'FontSize', 30, 'fontname','L M Roman10');
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
% figure;
% h = daboxplot(force,'xtlabels', {'Control','LocalSurround'},'whiskers',0,...
% 'scatter',1,'scattersize',250,'scatteralpha',0.6,'outliers',0,...
% 'colors',[.7 .7 .6],'linkline',1,'withinlines', 1);
% set(gca, 'FontSize', 30, 'fontname','L M Roman10');
% ylim([-0.6 .5])

figure;
plot(0:6, zeros(1, 7), 'k--','LineWidth',3);hold on
h = daviolinplot(data1_rel_mean,'colors',c,'box',3,...
'boxcolor','w','scatter',2,'jitter',0,'scattercolor','same',...
'smoothing', .1,'boxspacing',1.3, ...
'scattersize',400,'scatteralpha',0.7,...
'xtlabels', condition_names,'linkline',1,'withinlines', 1);
ylim([-0.6 .6])
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]);
set(gca, 'FontSize', 30, 'fontname','L M Roman10');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define table for statistical analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First color differences
for i=1:size(colors_rgb_ls,1)
    for j = 1:size(experiments_ill,1)
        lab1 = rgb2lab(experiments_ill(j, :).*colors_rgb_ls(i,:), ...
            'ColorSpace','linear-rgb','WhitePoint',[1 1 1]);
        lab2 = rgb2lab(colors_rgb_ls(i,:), 'ColorSpace','linear-rgb',...
            'WhitePoint',[1 1 1]);
        % lab1(1) = lab2(1);
        colorshift(i, j) = deltaE00(lab2', lab1');

        num_part = size(cci_surround{i}.CCI.ccip_neutral, 1);
        num_partc = size(cci_control{i}.CCI.ccip_neutral, 1);

        if isempty([cci_surround{i}.coloredMatch_reflected{:, j}])
            match = repmat([0 0 0]', 1, num_part);
        else
            match = ...
                reshape(([cci_surround{i}.coloredMatch_reflected{:, j}]), ...
                [], num_part);
        end
        
        if isempty([cci_control{i}.coloredMatch_reflected{:, j}])
            matchc = repmat([0 0 0]', 1, num_partc);
        else
            matchc = ...
                reshape(([cci_control{i}.coloredMatch_reflected{:, j}]), ...
                [], num_partc);
        end

        colorshift_control{i, j} = ...
            [deltaE00(match, repmat(lab1', 1, num_part))...
            deltaE00(matchc, repmat(lab1', 1, num_partc))];
        colorshift_localsurround{i, j} = ...
            [deltaE00(match, repmat(lab2', 1, num_part))...
            deltaE00(matchc, repmat(lab2', 1, num_partc))];

    end
end

for i = 1:length(cci_surround)

    switch cci_method
        case 1
            aux = cci_surround{i}.CCI.ccip(:, 1:end); ...
            aux2 = cci_control{i}.CCI.ccip(:, 1:end);
        case 2
            aux = cci_surround{i}.CCI.ccip_neutral(:, 1:end); ...
            aux2 = cci_control{i}.CCI.ccip_neutral(:, 1:end);
        case 3
            aux = cci_surround{i}.CCI.BRp(:, 1:end); ...
            aux2 = cci_control{i}.CCI.BRp(:, 1:end);
        case 4
            aux = cci_surround{i}.CCI.BRp_neutral(:, 1:end); ...
            aux2 = cci_control{i}.CCI.BRp_neutral(:, 1:end);
        case 5
            aux = cci_surround{i}.CCI.Karl(:, 1:end); ...
            aux2 = cci_control{i}.CCI.Karl(:, 1:end);
        case 6
            aux = cci_surround{i}.CCI.cci(:, 1:end); ...
            aux2 = cci_control{i}.CCI.cci(:, 1:end);
        case 7
            aux = cci_surround{i}.CCI.cci_neutral(:, 1:end); ...
            aux2 = cci_control{i}.CCI.cci_neutral(:, 1:end);
        case 8
            aux = cci_surround{i}.CCI.BR(:, 1:end); ...
            aux2 = cci_control{i}.CCI.BR(:, 1:end);
        case 9
            aux = cci_surround{i}.CCI.BR_neutral(:, 1:end); ...
            aux2 = cci_control{i}.CCI.BR_neutral(:, 1:end);
    end

    aux2p = [cci_surround{i}.Part_names;
            cci_control{i}.Part_names];
    
    aux_part = size(aux, 1) + size(aux2, 1);

    aux_part2 = size(aux2p, 1);

    aux_ill = [repmat("Neutral",size(aux, 1),1);...
        repmat("Blue",size(aux, 1),1);repmat("Yellow",size(aux, 1),1);...
        repmat("Green",size(aux, 1),1);repmat("Red",size(aux, 1),1);...
        repmat("Neutral",size(aux2, 1),1);...
        repmat("Blue",size(aux2, 1),1);repmat("Yellow",size(aux2, 1),1);...
        repmat("Green",size(aux2, 1),1);repmat("Red",size(aux2, 1),1)];
    
    local_control = [local_control;
        repmat("Suppression",length(aux(:)),1);
        repmat("Baseline",length(aux2(:)),1)];

    cci_all = [cci_all; aux(:); aux2(:)];

    diff = [diff;aux(:) - aux2(:);aux2(:) - aux2(:)];

    illuminants_all = [illuminants_all; aux_ill];

    experiment_all = [experiment_all; 
                      repmat(local{i},length(aux(:)),1);
                      repmat(local{i},length(aux2(:)),1)];

    aux_partcell = repmat(aux2p, 1,size(aux,2));
    participants_all = [participants_all; aux_partcell(:)];
    

    reflected = repmat(colors_rgb_ls(i, :), size(experiments_ill,1), 1) ...
        .* experiments_ill;
    
    reflected_lab = rgb2lab(reflected, 'ColorSpace','linear-rgb',...
        'WhitePoint',[1 1 1]);
    
    reflectance_lab = rgb2lab(colors_rgb_ls(i, :), 'ColorSpace',...
        'linear-rgb','WhitePoint',[1 1 1]);
    
    % colorshifts_control_all = [colorshifts_control_all; ...
    %     [colorshift_control{i, :}]'];

    colorshifts_localsurround_all = [colorshifts_localsurround_all; ...
        [colorshift_localsurround{i, :}]'];

    deltaerror = deltaE00(repmat(reflectance_lab, ...
        size(reflected_lab,1), 1)', reflected_lab');
    aux_shift = [repmat(deltaerror(1),size(aux,1),1);...
        repmat(deltaerror(2),size(aux,1),1);...
        repmat(deltaerror(3),size(aux,1),1);
        repmat(deltaerror(4),size(aux,1),1);
        repmat(deltaerror(5),size(aux,1),1);...
        repmat(deltaerror(1),size(aux2,1),1);...
        repmat(deltaerror(2),size(aux2,1),1);...
        repmat(deltaerror(3),size(aux2,1),1);
        repmat(deltaerror(4),size(aux2,1),1);
        repmat(deltaerror(5),size(aux2,1),1)];

    colorshifts_all = [colorshifts_all;aux_shift];

    [part, ill] = size(cci_surround{i}.targetCompetitorFit);
    aux_comp_pos = [];
    for k=1:ill
       aux_comp = ...
           cell2mat(vertcat(cci_surround{i}.targetCompetitorFit{:, k}));
       if k==1
           order_comp = -2:1:2;
       else
           order_comp = -3:1:1;
       end
       for kk=1:part
            aux_comp_pos = [aux_comp_pos;...
                min(interp1(aux_comp(kk,2:end),order_comp,...
                aux_comp(kk,1),'linear','extrap'), order_comp(end))];
       end
    end

    [part, ill] = size(cci_control{i}.targetCompetitorFit);
    % aux_comp_pos = [];
    for k=1:ill
       aux_comp = ...
           cell2mat(vertcat(cci_control{i}.targetCompetitorFit{:, k}));
       if k==1
           order_comp = -2:1:2;
       else
           order_comp = -3:1:1;
       end
       for kk=1:part
            aux_comp_pos = [aux_comp_pos;...
                min(interp1(aux_comp(kk,2:end),order_comp,...
                aux_comp(kk,1),'linear','extrap'), order_comp(end))];
       end
    end

    space_competitors = [space_competitors;aux_comp_pos];


    if outdoor
        scene = [scene; repmat("Outdoor",length(aux(:)) + length(aux2(:)), 1)];
    else
        scene = [scene; repmat("Indoor",length(aux(:))  + length(aux2(:)), 1)];
    end

    switch i
        case 1 % Khaki
            aux_neig = [repmat("None",size(aux, 1),1);...
                repmat("Opposite",size(aux, 1),1);...
                repmat("Neighbor",size(aux, 1),1);...
        repmat("Neighbor",size(aux, 1),1);repmat("Opposite",size(aux, 1),1);
        repmat("None",size(aux2, 1),1);...
        repmat("Opposite",size(aux2, 1),1);...
                repmat("Neighbor",size(aux2, 1),1);...
        repmat("Neighbor",size(aux2, 1),1);repmat("Opposite",size(aux2, 1),1);];

        case 2 % Rose
            aux_neig = [repmat("None",size(aux, 1),1);...
                repmat("Opposite",size(aux, 1),1);...
                repmat("Neighbor",size(aux, 1),1);...
        repmat("Opposite",size(aux, 1),1);repmat("Neighbor",size(aux, 1),1);
        repmat("None",size(aux2, 1),1);...
        repmat("Opposite",size(aux2, 1),1);...
                repmat("Neighbor",size(aux2, 1),1);...
        repmat("Opposite",size(aux2, 1),1);repmat("Neighbor",size(aux2, 1),1);];

        case 3 % Purple
            aux_neig = [repmat("None",size(aux, 1),1);...
                repmat("Neighbor",size(aux, 1),1);...
                repmat("Opposite",size(aux, 1),1);...
        repmat("Opposite",size(aux, 1),1);repmat("Neighbor",size(aux, 1),1);
        repmat("None",size(aux2, 1),1);...
        repmat("Neighbor",size(aux2, 1),1);...
                repmat("Opposite",size(aux2, 1),1);...
        repmat("Opposite",size(aux2, 1),1);repmat("Neighbor",size(aux2, 1),1);];

        case 4 % Teal
            aux_neig = [repmat("None",size(aux, 1),1);...
                repmat("Neighbor",size(aux, 1),1);...
                repmat("Opposite",size(aux, 1),1);...
        repmat("Neighbor",size(aux, 1),1);repmat("Opposite",size(aux, 1),1);
        repmat("None",size(aux2, 1),1);...
        repmat("Neighbor",size(aux2, 1),1);...
                repmat("Opposite",size(aux2, 1),1);...
        repmat("Neighbor",size(aux2, 1),1);repmat("Opposite",size(aux2, 1),1);];

    end
    neighbor_all = [neighbor_all; aux_neig(:)];
end

tablenames = ...
    {'CCI', 'CCIdiff', 'MatchPosition', 'Participant', 'LocalSurround','Baseline', ...
    'ColorShift', 'ColorShiftLS', ...
    'Illuminant', 'Direction', 'Scene'};
tbl = table(cci_all, diff, space_competitors,...
    participants_all, experiment_all, local_control, ...
    colorshifts_localsurround_all, colorshifts_all, ...
    illuminants_all, neighbor_all, scene, 'VariableNames', tablenames);

tbl.LocalSurround = nominal(tbl.LocalSurround);
tbl.Illuminant    = nominal(tbl.Illuminant);
tbl.Direction     = nominal(tbl.Direction);
tbl.Baseline      = nominal(tbl.Baseline);
tbl.Scene         = nominal(tbl.Scene);

tbl_neutral = tbl(tbl.Illuminant == 'Neutral', :);
tbl_colored = tbl;
tbl_colored(tbl.Illuminant == 'Neutral', :) = [];

writetable(tbl, strcat(save_table_filename, '.csv'));
writetable(tbl_neutral, strcat(save_table_filename, '_neutral.csv'));
writetable(tbl_colored, strcat(save_table_filename, '_colored.csv'));

%% Only the colour illuminants
% disp('Only coloured illuminants')
% tbl = tbl_colored;

%% Statistical analysis based on the direction: Neighboring, Opposing
lme{1} = fitlme(tbl, 'CCI ~ 1 + Direction + (1|Participant)');  
lme{2} = fitlme(tbl, 'CCI ~ 1 + Baseline + (1|Participant)');  
lme{3} = fitlme(tbl, ['CCI ~ 1 + ' ...
    'LocalSurround*Illuminant + (1|Participant)']); 
% lme{4} = fitlme(tbl, ['CCI ~ 1 + ' ...
%     'LocalSurround*Illuminant + Baseline + (1|Participant)']); 

glme{1} = fitglme(tbl, 'CCI ~ 1 + Direction + (1|Participant)');
glme{2} = fitglme(tbl, 'CCI ~ 1 + Baseline + (1|Participant)'); 
glme{3} = fitglme(tbl, ['CCI ~ 1 + ' ...
    'LocalSurround*Illuminant + (1|Participant)']); 

lme_dir = anova(lme{1});
lme_base = anova(lme{2});

glme_dir = anova(glme{1});
glme_base = anova(glme{2});

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
lme_ls{1} = fitlme(tbl_colored, ['CCIdiff ~ 1 + LocalSurround + ' ...
    '(1|Participant)']);  
lme_ls{2} = fitlme(tbl_colored, ['CCI ~ 1 + LocalSurround + ' ...
    '(1|Participant)']); 
 
ls_ccidiff = anova(lme_ls{1});
ls_cci = anova(lme_ls{2});
disp(ls_cci)
disp(ls_ccidiff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot CCI results based on the neighboring, opposing illuminants
c2 = [125,125,155;
     176,196,222;
     112,128,144;
     125,125,125;
     176,196,222]./255;

data2{1} = tbl.CCI(tbl.Direction == 'Opposite' & ...
    tbl.Baseline == 'Baseline');
data2{2} = tbl.CCI(tbl.Direction == 'Neighbor' & ...
    tbl.Baseline == 'Baseline');
data2{3} = tbl.CCI(tbl.Direction == 'Opposite' & ...
    tbl.Baseline == 'Suppression');
data2{4} = tbl.CCI(tbl.Direction == 'Neighbor' & ...
    tbl.Baseline == 'Suppression');

if relative
    data3{1} = -1 + (data2{1}) ./ data2{1};
    data3{2} = -1 + (data2{2}) ./ data2{2};
    data3{3} = -1 + (data2{3}) ./ data2{1}; 
    data3{4} = -1 + (data2{4}) ./ data2{2};
else
    data3{1} = data2{1} - data2{1};
    data3{2} = data2{2} - data2{2};
    data3{3} = data2{3} - data2{1};
    data3{4} = data2{4} - data2{2};
end

unique_part = convertCharsToStrings(unique(participants_all));

for i=1:length(unique_part)
    
    opp(i) = mean(tbl.CCI(tbl.Direction == 'Opposite' & ...
        tbl.Baseline == 'Suppression' & ...
        convertCharsToStrings(tbl.Participant) == unique_part(i, :)));
    neig(i) = mean(tbl.CCI(tbl.Direction == 'Neighbor' & ...
        tbl.Baseline == 'Suppression' & ...
        convertCharsToStrings(tbl.Participant) == unique_part(i, :)));

    dnoneo(i) = mean(tbl.CCI(tbl.Direction == 'Opposite' & ...
        tbl.Baseline == 'Baseline' & ...
        convertCharsToStrings(tbl.Participant) == unique_part(i, :)));

    dnonen(i) = mean(tbl.CCI(tbl.Direction == 'Neighbor' & ...
        tbl.Baseline == 'Baseline' & ...
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
            if ~isempty(tbl.CCIdiff(tbl.Baseline == 'Suppression' & ...
                    tbl.LocalSurround == ls & ...
                    convertCharsToStrings(tbl.Participant) == unique_part(i, :) & ...
                    tbl.Illuminant == illuminants_names(k)))

                all_ls(i, j, k) = tbl.CCIdiff(tbl.Baseline == 'Suppression' & ...
                    tbl.LocalSurround == ls & ...
                    convertCharsToStrings(tbl.Participant) == unique_part(i, :) & ...
                    tbl.Illuminant == illuminants_names(k));
            end

            if ~isempty(tbl.CCI(tbl.Baseline == 'Baseline' & ...
                    tbl.LocalSurround == ls & ...
                    convertCharsToStrings(tbl.Participant) == unique_part(i, :) & ...
                    tbl.Illuminant == illuminants_names(k)))

                all_c_ls(i, j, k) = tbl.CCI(tbl.Baseline == 'Baseline' & ...
                    tbl.LocalSurround == ls & ...
                    convertCharsToStrings(tbl.Participant) == unique_part(i, :) & ...
                    tbl.Illuminant == illuminants_names(k));
            end

        end
        
    end

end
parti_local = nanmean(all_ls, 3);
parti_control = nanmean(all_c_ls, 3);
for i=1:length(unique_part)
    for j=i+1:length(unique_part)
        aux = corrcoef(parti_local(i, :), parti_local(j, :), ...
            'rows', 'complete');
        Cinter(i, j) = aux(1, 2);
        Cinter(j, i) = aux(1, 2);

        aux = corrcoef(parti_control(i, :), parti_control(j, :), ...
            'rows', 'complete');
        Cinter_c(i, j) = aux(1, 2);
        Cinter_c(j, i) = aux(1, 2);
    end
end

parti = nanmean(nanmean(all_ls, 3),2);
parti_c = nanmean(nanmean(all_c_ls, 3),2);
figure;
plot(0:2, zeros(3, 1), '--k', 'LineWidth', 3);hold on
h = daboxplot(parti,...
    'xtlabels', {'Suppression'},'whiskers',0,...
'scatter',1,'scattersize',350,'scatteralpha',0.6,'withinlines',1,...
'linkline',1,'outliers',0,'colors',c);
set(gca, 'FontSize', 30, 'fontname','L M Roman10');
axis([0 2 -.2 .1])
% h = daboxplot([parti parti_c],...
%     'xtlabels', {'Suppression','Baseline'},'whiskers',0,...
% 'scatter',1,'scattersize',350,'scatteralpha',0.6,'withinlines',1,...
% 'linkline',1,'outliers',0,'colors',c);
% set(gca, 'FontSize', 30, 'fontname','L M Roman10');

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
h = daboxplot(data7,'xtlabels', {''},'whiskers',0,...
'scatter',1,'scattersize',350,'scatteralpha',0.6,'withinlines',1,...
'linkline',1,'outliers',0,'colors',c2);
set(gca, 'FontSize', 25, 'fontname','L M Roman10');
ylim([-.5 .5])
ylabel('Absolute error (\Delta AE)')

% figure;h2 = daviolinplot(data2,'colors',c2,'box',1,...
% 'boxcolor','w','scatter',2,'jitter',1,'scattercolor','same',...
% 'smoothing', .15,'boxspacing',1.3, ...
% 'scattersize',200,'scatteralpha',0.7,'linkline',1,...
% 'xtlabels', {"ControlO","ControlN","Opposing","Neighboring"}); % condition_names2

