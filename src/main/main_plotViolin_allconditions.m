clear
% close all

reflected = 1;
outdoor = 1;
median_sel = 0;
relative = 1;

N = (2*reflected+outdoor)^(2*reflected+1) + 3*median_sel;

colors_surround = {'Khaki', 'Rose', 'Purple', 'Teal'};

colors_ill_clockwise = {'Yellow', 'Red','Blue', 'Green'};
colors_ill_anticlockwise = {'Green', 'Yellow', 'Red','Blue'};
colors_ill_oppclockwise = {'Blue', 'Green','Yellow', 'Red'};
colors_ill_oppanticlockwise = {'Red', 'Blue', 'Green','Yellow'};

% order_illuminants = {'Blue', 'Green', 'Neutral', 'Yellow', 'Red'};
order_illuminants = {'Neutral','Blue', 'Yellow','Green', 'Red'};

plot_labels = {'Blue', 'Green','Yellow','Red', ...
    'Blue', 'Green','Yellow','Red', ...
    'Blue', 'Green','Yellow','Red', ...
    'Blue', 'Green','Yellow','Red'};


switch(N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MLDS selection for the MATCH    
    case 0 % 0 reflected 0 outdoor
        cci_surround{1} = load('Kaki_indoor_reflectance11_cci.mat');
        cci_surround{2} = load('Rose_indoor_reflectance11_cci.mat');
        cci_surround{3} = load('Purple_indoor_reflectance11_cci.mat');
        cci_surround{4} = load('Teal_indoor_reflectance11_cci.mat');
        cci_surround_control = load('Control_indoor_reflectance11_cci.mat');
    case 8 % 1 reflected 0 outdoor
        cci_surround{1} = load('Kaki_indoor_reflectance11_cci.mat');
        cci_surround{2} = load('Rose_indoor_reflectance11_cci.mat');
        cci_surround{3} = load('Purple_indoor_reflectance11_cci.mat');
        cci_surround{4} = load('Teal_indoor_reflectance11_cci.mat');
        cci_surround_control = load('Control_indoor_reflectance11_cci.mat');
    case 1 % 0 reflected 1 outdoor
        cci_surround{1} = load('Kaki_outdoor_reflectance11_cci.mat');
        cci_surround{2} = load('Rose_outdoor_reflectance11_cci.mat');
        cci_surround{3} = load('Purple_outdoor_reflectance11_cci.mat');
        cci_surround{4} = load('Teal_outdoor_reflectance11_cci.mat');
        cci_surround_control = load('Control_outdoor_reflectance11_cci.mat');
    case 27 % 1 reflected 1 outdoor
        cci_surround{1} = load('khaki_outdoor.mat');
        cci_surround{2} = load('rose_outdoor.mat');
        cci_surround{3} = load('purple_outdoor.mat');
        cci_surround{4} = load('teal_outdoor.mat');
        cci_surround_control = load('control_khaki_outdoor.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Median selection for the MATCH

   case 3 % 0 reflected 0 outdoor
        cci_surround{1} = load('Kaki_indoor_median_cci.mat');
        cci_surround{2} = load('Rose_indoor_median_cci.mat');
        cci_surround{3} = load('Purple_indoor_median_cci.mat');
        cci_surround{4} = load('Teal_indoor_median_cci.mat');
        cci_surround_control = load('Control_indoor_median_cci.mat');
    case 11 % 1 reflected 0 outdoor
        cci_surround{1} = load('Kaki_indoor_median_reflected_cci.mat');
        cci_surround{2} = load('Rose_indoor_median_reflected_cci.mat');
        cci_surround{3} = load('Purple_indoor_median_reflected_cci.mat');
        cci_surround{4} = load('Teal_indoor_median_reflected_cci.mat');
        cci_surround_control = load('Control_indoor_median_reflected_cci.mat');
    case 4 % 0 reflected 1 outdoor
        cci_surround{1} = load('Kaki_outdoor_median_cci.mat');
        cci_surround{2} = load('Rose_outdoor_median_cci.mat');
        cci_surround{3} = load('Purple_outdoor_median_cci.mat');
        cci_surround{4} = load('Teal_outdoor_median_cci.mat');
        cci_surround_control = load('Control_outdoor_median_cci.mat');
    case 30 % 1 reflected 1 outdoor
        cci_surround{1} = load('Kaki_outdoor_median_reflected_cci.mat');
        cci_surround{2} = load('Rose_outdoor_median_reflected_cci.mat');
        cci_surround{3} = load('Purple_outdoor_median_reflected_cci.mat');
        cci_surround{4} = load('Teal_outdoor_median_reflected_cci.mat');
        cci_surround_control = load('Control_outdoor_median_reflected_cci.mat');
end


mean_surround_control_ill = mean(...
    cci_surround_control.CCI.cci_neutral(:, 2:end));                           % cci_surround_control.control.pert_cc([1:2 4:5], :));

cci_all = [];
participants_all = [];
experiment_all = [];
illuminants_all = [];
illuminants_all_rgb = [];
colorshifts_all = [];

aux_ill = [repmat("blue",1,1);repmat("green",1,1);...
        repmat("yellow",1,1);repmat("red",1,1)];
local = {"Khaki", "Rose","Purple","Teal","Control"};
localcolor = [0.2929    0.2796    0.1624;0.4782    0.2157    0.2201;...
    0.3055    0.2547    0.3785;0.1989    0.2962    0.2960;.2 .2 .2];
experiments_ill = [0.2316 0.2833 0.4579; %blue
    0.233 0.3 0.2003;  %green
    0.38 0.2441 0.1305;  %yellow
    0.4152 0.2114 0.3961]; %red
%    0.2912 0.2705 0.2614; %neutral

c = [125,125,125;
    189,183,107;
    205,92,92;
    147,112,219;
    95,158,160]./255;

data1{1} = mean(cci_surround_control.CCI.cci_neutral(:, 2:end)', 2);
data1{2} = mean(cci_surround{1}.CCI.cci_neutral(:, 2:end)', 2);
data1{3} = mean(cci_surround{2}.CCI.cci_neutral(:, 2:end)', 2);
data1{4} = mean(cci_surround{3}.CCI.cci_neutral(:, 2:end)', 2);
data1{5} = mean(cci_surround{4}.CCI.cci_neutral(:, 2:end)', 2);

condition_names = {'Water', 'Land', 'Moon', 'Hyperspace'};

figure;h = daviolinplot(data1,'colors',c,'box',3,...
'boxcolor','w','scatter',2,'jitter',0,'scattercolor','same',...
'smoothing', .1,'boxspacing',1.3, ...
'scattersize',120,'scatteralpha',0.7,...
'xtlabels', condition_names);
ylim([0 1.2])
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]);

for i=1:size(localcolor,1)
    for j = 1:size(experiments_ill,1)
        lab1 = rgb2lab(experiments_ill(j, :).*localcolor(i,:), ...
            'ColorSpace','linear-rgb','WhitePoint',[1 1 1]);
        lab2 = rgb2lab(localcolor(i,:), 'ColorSpace','linear-rgb',...
            'WhitePoint',[1 1 1]);
        % lab2(1) = lab1(1);
        colorshift(i, j) = deltaE00(lab2', lab1');
    end
end
%colorshift(i+1, :) = ones(1,size(colorshift,2));
colorshift = colorshift';
%% Add all the conditions together
aux_control = cci_surround_control.CCI.cci_neutral(:, 2:end);
for i=1:length(cci_surround)

    % if relative && ~reflected
    %     aux = 100 - (aux_control - cci_surround{i}.CCI.cci_neutral(:, 2:end));
    % elseif relative && reflected
    %     aux = 1 - (aux_control - cci_surround{i}.CCI.cci_neutral(:, 2:end));
    % else
    %     aux = cci_surround{i}.pert_cc([1:2 4:5], :);
    % end
    aux = cci_surround{i}.CCI.cci_neutral(:, 2:end)';
    cci_all = [cci_all; ...
        aux(:)];

    for j=1:size(aux,2)

        if j <= 10
            participants_all = [participants_all;...
                repmat(strcat("p00", num2str(j-1)),4,1)];
        else
            participants_all = [participants_all;...
                repmat(strcat("p0", num2str(j-1)),4,1)];
        end
        illuminants_all = [illuminants_all;aux_ill];
        illuminants_all_rgb = [illuminants_all_rgb;experiments_ill];
        colorshifts_all = [colorshifts_all;colorshift(:, i)];
    end
    experiment_all = [experiment_all;repmat(local{i},length(aux(:)),1)];
    
end


tbl = table(cci_all,experiment_all, participants_all,illuminants_all,colorshifts_all);
lme = fitlme(tbl,'cci_all ~ illuminants_all + (illuminants_all|participants_all)');
lme2 = fitlme(tbl,'cci_all ~ experiment_all + illuminants_all + (experiment_all|participants_all) + (illuminants_all|participants_all)');

% figure;h = daviolinplot(1-(data1{1} - data1),'colors',c,'box',3,...
% 'boxcolor','w','scatter',2,'jitter',0,'scattercolor','same',...
% 'smoothing', .1,'boxspacing',1.3, ...
% 'scattersize',120,'scatteralpha',0.7,...
% 'xtlabels', condition_names);
% ylim([0 1.2])
% xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]);

if relative && ~reflected
    aux2 = 100 - (aux_control - aux_control);
elseif relative && reflected
    aux2 = 1 - (aux_control - aux_control);
else
    aux2 = aux_control;
end
cci_all = [cci_all; ...
    aux2(:)];
for j=1:size(aux,2)
    if j <= 10
        participants_all = [participants_all;...
            repmat(strcat("p00", num2str(j-1)),4,1)];
    else
        participants_all = [participants_all;...
            repmat(strcat("p0", num2str(j-1)),4,1)];
    end

    illuminants_all = [illuminants_all;aux_ill];
    illuminants_all_rgb = [illuminants_all_rgb;experiments_ill];
    colorshifts_all = [colorshifts_all;colorshift(:, i+1)];

end
experiment_all = [experiment_all;repmat(local{i+1},length(aux(:)),1)];

factors_all = {experiment_all, illuminants_all, ...
    participants_all, colorshifts_all};

aov = anova(factors_all,cci_all(:),CategoricalFactors=[1 2 3],...
    RandomFactors=3,...
    FactorNames=["Experiments","Illuminants",...
    "Participants","ColorShift"]) % ,ModelSpecification="interactions"

p = anovan(cci_all(:),{experiment_all illuminants_all ...
    participants_all colorshifts_all},...
    'model','linear',...
    'random', 3, ...
    'continuous', 4, ...
    'varnames',{'experiment','illuminants','participants', 'colorshift'})

% p = anovan([cci_all(:);cci_all(:);cci_all(:);cci_all(:)],...
%     {[experiment_all;experiment_all;experiment_all;experiment_all]...
%     [illuminants_all;illuminants_all;illuminants_all;illuminants_all] ...
%     [participants_all;participants_all;participants_all;participants_all]...
%     [colorshifts_all;colorshifts_all;colorshifts_all;colorshifts_all]},...
%     'model','full',...
%     'varnames',{'experiment','illuminants','participants', 'colorshift'})

[pvals,tbl,stats] = anovan(cci_all(:), {illuminants_all participants_all}, ... 
'model',2, 'random',2,'varnames',{'illuminants','participants'});

figure;
tiledlayout(2,2)
nexttile
boxchart(aov,["Experiments","Participants"])
ylabel("CCI")
xlabel("Experiments")
nexttile
boxchart(aov,["Experiments","Illuminants"])
ylabel("CCI")
xlabel("Experiments")
nexttile
boxchart(aov,"Participants")
ylabel("CCI")
xlabel("Participants")

figure;
tiledlayout(2,2)
nexttile
boxchart(aov,["Experiments"])
ylabel("CCI")
xlabel("Experiments")
nexttile
boxchart(aov,["Illuminants"])
ylabel("CCI")
xlabel("Illuminants")
nexttile
boxchart(aov,"Participants")
ylabel("CCI")
xlabel("Participants")

figure;
%%% Control and surrounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(local)
    disp(local{i})
    data = cci_all(strcmp(experiment_all, local{i}) );
    factor_ill = illuminants_all( strcmp(experiment_all, local{i}) );
    factor_part = participants_all( strcmp(experiment_all, local{i}) );
    factor_exp = experiment_all( strcmp(experiment_all, local{i}) );
    factor_shift = colorshifts_all( strcmp(experiment_all, local{i}) );
    factor_ill_rgb = ...
        illuminants_all_rgb( strcmp(experiment_all, local{i}), : );

    n = grp2idx(factor_ill);

    anova(factor_ill,data(:),CategoricalFactors=1,...
        FactorNames=["Illuminants"]) % ,ModelSpecification="interactions"

    pp(i) = anovan(data(:),{factor_ill},...
        'model','linear',...
        'varnames',{'illuminants'})

    % factor_ill = nominal(factor_ill);
    % factor_part = nominal(factor_part);
    % factor_exp = nominal(factor_exp);

    tbl = table(data,factor_ill,factor_part);
    %Fit a linear regression model 

    % mdl = fitlm(tbl,'data ~ 0+factor_ill')
    lme = fitlme(tbl,'data ~ 1 + factor_ill + (1|factor_part)')

    tbl = table(data,factor_ill,factor_part,factor_exp,factor_shift);
    % lme2 = fitlme(tbl,['data ~ factor_shift +  (1|factor_part)'])
    % lme = fitlme(tbl,...
    %     ['data ~ factor_ill*factor_exp + (1|factor_part)'])
    

    mdl = fitglme(tbl, 'data ~ factor_ill + (1|factor_part)');
    emm = emmeans(mdl, {'factor_ill'}, 'unbalanced');
    % (after inspecting emm.table and selecting your contrast)
    pVal = contrasts_wald(mdl, emm, [1 -1 0 0]);                            %% other combinations [1 0 -1 0]

    % scatter(factor_shift, data, 100,'filled', 'MarkerFaceColor',...
    %     localcolor(i, :).^(.45)); hold on; %  factor_ill_rgb
    auxdata = [mean(data(n==1)) mean(data(n==2)) ...
        mean(data(n==3)) mean(data(n==4))];
    auxshift = [mean(factor_shift(n==1)) mean(factor_shift(n==2)) ...
        mean(factor_shift(n==3)) mean(factor_shift(n==4))];
    % scatter3(auxshift, unique(n), auxdata,...
    %     150, 'filled', 'MarkerFaceColor',...
    %     localcolor(i, :).^(.4)); hold on; %  factor_ill_rgb
    
    scatter3(factor_shift, n, data,...
        data.*200, 'filled', 'MarkerFaceColor',...
        localcolor(i, :).^(.4)); hold on; %  factor_ill_rgb
    
    % subplot(1,3,1);
    % scatter(factor_shift, n,...
    %     data.*100, 'filled', 'MarkerFaceColor',...
    %     localcolor(i, :).^(.4)); hold on; %  factor_ill_rgb
    % subplot(1,3,2);
    % scatter(factor_shift, data,...
    %     n.*100, 'filled', 'MarkerFaceColor',...
    %     localcolor(i, :).^(.4)); hold on; %  factor_ill_rgb
    % subplot(1,3,3);
    % scatter(n, data,...
    %     factor_shift.*10, 'filled', 'MarkerFaceColor',...
    %     localcolor(i, :).^(.4)); hold on; %  factor_ill_rgb

    unique_part = unique(factor_part);
    for j=1:length(unique_part)
        data_aux(j) = mean(tbl.data(tbl.factor_part == unique_part{j}, :));
        part_aux{j} = unique_part{j};
        exp_aux{j} = factor_exp{1};
    end
    tbl_local{i} = table(data_aux', part_aux', exp_aux', ...
        'VariableNames', {'data','participants','experiment'});
end
yticks([1 2 3 4])
yticklabels({'Blue', 'Green','Yellow', 'Red'})
xlabel('Shift in DeltaE')
ylabel('Illuminant colors')
zlabel('Color Constancy Index')
set(gca, 'FontSize', 20, 'fontname','Times New Roman');


labels = {'Shift' 'Illuminant' 'CCI'};
d = [factor_shift, n, data];

[h,ax] = plotmatrix(d);                        % create a 4 x 4 matrix of plots
for i = 1:3                                      % label the plots
  xlabel(ax(3,i), labels{i})
  ylabel(ax(i,1), labels{i})
end


cc = [33 33 33;0, 131, 143;106, 27, 154;229, 115, 115;158, 157, 36]./255;

figure;
if relative
    plot(0.5:.5:4.5, ones(1, length(0.5:.5:4.5)), 'k--', 'LineWidth',2.5);...
        hold on
end
if relative
    start = length(local)-1;
    count = 2;
else
    start = length(local);
    count = 1;
end
exp_num = 1;
for i = start:-1:1
    
    if ~reflected
        tbl_local{i}.data = tbl_local{i}.data./100;
    end
    vs = Violin({tbl_local{i}.data},...
        exp_num,...
        'HalfViolin','full',...% left, full
        'QuartileStyle','shadow',... % boxplot, none
        'DataStyle', 'scatter',... % scatter, none
        'ShowNotches', false,...
        'ShowMean', true,...
        'ShowMedian', true,...
        'MedianMarkerSize', 96,...
        'MarkerSize',80, ...
        'LineWidth', 3,...
        'ViolinColor', {cc(count, :)});
    rev_local{exp_num} = local{i};
    count = count +1;
    exp_num = exp_num + 1;
    
end

if relative
    xticks([1 2 3 4])
    axis([0.5 4.5 0.5 1.3])
    ylabel('CCI relative to Control')
else
    xticks([1 2 3 4 5])
    axis([0.5 5.5 0 1.3])
    ylabel('Color Constancy Index')

end
xticklabels(rev_local)

set(gca, 'FontSize', 34, 'fontname','Times New Roman');

[~,p_k] = ttest(tbl_local{5}.data, ...
    tbl_local{1}.data);
[~,p_r] = ttest(tbl_local{5}.data, ...
    tbl_local{2}.data);
[~,p_t] = ttest(tbl_local{5}.data, ...
    tbl_local{4}.data);
[~,p_p] = ttest(tbl_local{5}.data, ...
    tbl_local{3}.data);


cci_all = cci_all(end:-1:1);
illuminants_all = illuminants_all(end:-1:1);
participants_all = participants_all(end:-1:1);...
experiment_all = experiment_all(end:-1:1);
colorshifts_all = colorshifts_all(end:-1:1);


tbl = table(cci_all,illuminants_all,participants_all,...
    experiment_all,colorshifts_all);


mdl2 = fitglme(tbl, ...
    'cci_all ~ experiment_all + (1|participants_all)');
aa = anova(mdl2);
emm2 = emmeans(mdl2, {'experiment_all'}, 'unbalanced');
pVal = contrasts_wald(mdl2, emm2, [-1 -1 -1 -1 1]);

[~,p_k] = ttest(tbl.cci_all(tbl.experiment_all == 'Control'), ...
    tbl.cci_all(tbl.experiment_all == 'Khaki'));
[~,p_r] = ttest(tbl.cci_all(tbl.experiment_all == 'Control'), ...
    tbl.cci_all(tbl.experiment_all == 'Rose'));
[~,p_t] = ttest(tbl.cci_all(tbl.experiment_all == 'Control'), ...
    tbl.cci_all(tbl.experiment_all == 'Teal'));
[~,p_p] = ttest(tbl.cci_all(tbl.experiment_all == 'Control'), ...
    tbl.cci_all(tbl.experiment_all == 'Purple'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Individual differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

part = unique(participants_all);
corrmat = nan(length(part));
for i=1:length(part)
    aux1 = aov.Y(strcmp(aov.Factors.Participants, part{i}));
    for j=i+1:length(part)
        aux2 = aov.Y(strcmp(aov.Factors.Participants, part{j}));
        aux3 = corrcoef(aux1,aux2);
        corrmat(i,j) = aux3(1,2);
    end

end
figure;
heatmap(corrmat);
colormap(parula)
disp('Mean Individual differences')
mean(mean(corrmat', 'omitnan'),'omitnan')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=length(local):-1:1
    Y1{:, length(local) - i + 1} = aov.Y(strcmp(aov.Factors.Experiments,...
        local{i}) );
    xticks_labels{length(local) - i + 1} = local{i};
end
if ~reflected
    for i=length(local):-1:1
        Y1{:, length(local) - i + 1} = aov.Y(...
            strcmp(aov.Factors.Experiments, local{i}) )./100;
    end
end

cc = [33 33 33;0, 131, 143;106, 27, 154;229, 115, 115;158, 157, 36]./255;

figure;
for i = 1:length(local)
    vs = Violin({Y1{:, i}},...
        i,...
        'HalfViolin','full',...% left, full
        'QuartileStyle','shadow',... % boxplot, none
        'DataStyle', 'scatter',... % scatter, none
        'ShowNotches', false,...
        'ShowMean', true,...
        'ShowMedian', true,...
        'MedianMarkerSize', 96,...
        'MarkerSize',80, ...
        'LineWidth', 3,...
        'ViolinColor', {cc(i, :)});
end
axis([0.5 5.5 0 1.3])
xticks([1 2 3 4 5])
xticklabels(xticks_labels)
ylabel('Color Constancy Index')
set(gca, 'FontSize', 44, 'fontname','Times New Roman');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Neighbors clockwise
order_ill_clockwise = ...
    names2numbers(order_illuminants, colors_ill_clockwise);
for i=1:length(cci_surround)
    cci_surround_clockwise(:, order_ill_clockwise(i)) = ...
        cci_surround{i}.CCI.cci_neutral(:, order_ill_clockwise(i));
end

%% Neighbors anticlockwise
order_ill_anticlockwise = ...
    names2numbers(order_illuminants, colors_ill_anticlockwise);
for i=1:length(cci_surround)
    cci_surround_anticlockwise(:, order_ill_anticlockwise(i)) = ...
        cci_surround{i}.CCI.cci_neutral(:, order_ill_clockwise(i));
end

%% Opposite clockwise
order_ill_oppclockwise = ...
    names2numbers(order_illuminants, colors_ill_oppclockwise);
for i=1:length(cci_surround)
    cci_surround_oppclockwise(:, order_ill_oppclockwise(i)) = ...
        cci_surround{i}.CCI.cci_neutral(:, order_ill_clockwise(i));
end

%% Opposite anticlockwise
order_ill_oppanticlockwise = ...
    names2numbers(order_illuminants, colors_ill_oppanticlockwise);
for i=1:length(cci_surround)
    cci_surround_oppanticlockwise(:, order_ill_oppanticlockwise(i)) = ...
        cci_surround{i}.CCI.cci_neutral(:, order_ill_clockwise(i));
end

%% Plot all illuminants
ci_colorimeter_direct = [cci_surround_clockwise([1:2 4:5], :)'...
    cci_surround_anticlockwise([1:2 4:5], :)' ... 
    cci_surround_oppclockwise([1:2 4:5], :)' ...
    cci_surround_oppanticlockwise([1:2 4:5], :)'];

ci_colorimeter_direct_c = [cci_surround_control.pert_cc([1:2 4:5], :)'...
    cci_surround_control.pert_cc([1:2 4:5], :)' ... 
    cci_surround_control.pert_cc([1:2 4:5], :)' ...
    cci_surround_control.pert_cc([1:2 4:5], :)'];

aa = ci_colorimeter_direct';

if relative
    aa = 1 - (ci_colorimeter_direct_c' - ci_colorimeter_direct');
end

participants = [];
experiment = [];
illuminants = [];

for i=1:size(aa,2)
    participants = [participants;repmat(strcat("p0", num2str(i-1)),16,1)];
    experiment = [experiment;repmat("neighboring",4,1);...
        repmat("neighboring",4,1);repmat("opposite",4,1);...
        repmat("opposite",4,1)];
    illuminants = [illuminants;repmat("blue",1,1);repmat("green",1,1);...
        repmat("yellow",1,1);repmat("red",1,1);repmat("blue",1,1);...
        repmat("green",1,1);repmat("yellow",1,1);repmat("red",1,1);...
        repmat("blue",1,1);repmat("green",1,1);...
        repmat("yellow",1,1);repmat("red",1,1);repmat("blue",1,1);...
        repmat("green",1,1);repmat("yellow",1,1);repmat("red",1,1)];
end


factors = {experiment, illuminants, participants}
aov = anova(factors,aa(:),FactorNames=["Experiments","Illuminants",...
"Participants"])
Y3{:, 1} = aov.Y(strcmp(aov.Factors.Experiments, 'neighboring') );
Y3{:, 2} = aov.Y(strcmp(aov.Factors.Experiments, 'opposite') );
if ~reflected
    Y3{:, 1} = aov.Y(strcmp(aov.Factors.Experiments, 'neighboring') )./100;
    Y3{:, 2} = aov.Y(strcmp(aov.Factors.Experiments, 'opposite') )./100;
end

[H,P] = ttest(Y3{1,1},Y3{1,2})
part_unique = unique(aov.Factors.Participants);

for i = 1:length(part_unique)

    Y4(i, 1) = mean(aov.Y(aov.Factors.Participants == part_unique{i} & ...
        aov.Factors.Experiments == 'neighboring'));
    Y4(i, 2) = mean(aov.Y(aov.Factors.Participants == part_unique{i} & ...
        aov.Factors.Experiments == 'opposite'));
    if ~reflected
        Y4(i, 1) = Y4(i, 1) ./ 100;
        Y4(i, 2) = Y4(i, 2) ./ 100;
    end
end

% figure;
% violin(Y3,'x', [1 3], 'facecolor',[.3 .3 .3;0.9765 0.6588 0.1451],'edgecolor',...
%     'none','bw',0.1,'mc','k','medc','r-.')
% axis([0 4 0 1.3])
% xticks([1 3])
% xticklabels({'Neighboring','Opposite'})
% set(gca, 'FontSize', 24, 'fontname','Times New Roman');
% ylabel('Color Constancy Index')

cc1 = [.3 .3 .3;0.9765 0.6588 0.1451];
cc1 = [69 90 100;96 125 139]./255;
figure;
if relative
    plot(0.5:2.5, ones(length(0.5:2.5), 1), 'k--','LineWidth',2.5);hold on
end
for i = 1:length(Y3)
    vs = Violin({Y3{:, i}},...
        i,...
        'HalfViolin','full',...% left, full
        'QuartileStyle','shadow',... % boxplot, none
        'DataStyle', 'scatter',... % scatter, none
        'ShowNotches', false,...
        'ShowMean', true,...
        'ShowMedian', true,...
        'MedianMarkerSize', 96,...
        'MarkerSize',80, ...
        'ViolinAlpha', {.3, .9},...
        'LineWidth', 3, ...
        'ViolinColor', {cc1(i, :)});
end
axis([0.5 2.5 0.5 1.3])
xticks([1 2])
xticklabels({'Neighboring illuminants','Opposite illuminants'})
if relative
    ylabel('CCI relative to Control')
else
    ylabel('Color Constancy Index')

end
set(gca, 'FontSize', 24, 'fontname','Times New Roman');

c = [125, 125, 125;
    70,130,180;
95,158,160]./255;

data3{1} = data1{1};
data3{2} = Y4(:, 1);
data3{3} = Y4(:, 2);

figure;h = daviolinplot(Y3,'colors',c,'box',3,...
'boxcolor','w','scatter',2,'jitter',0,'scattercolor','same',...
'smoothing', .1,'boxspacing',1.3, ...
'scattersize',120,'scatteralpha',0.7,...
'xtlabels', condition_names);
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]);

figure;h = daviolinplot(Y4,'colors',c,'box',3,...
'boxcolor','w','scatter',2,'jitter',0,'scattercolor','same',...
'smoothing', .1,'boxspacing',1.3, ...
'scattersize',150,'scatteralpha',0.7,...
'xtlabels', condition_names);
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]);
ylim([0 1.4])

figure;
if relative
    plot(0.5:2.5, ones(length(0.5:2.5), 1), 'k--','LineWidth',2.5);hold on
end
for i = 1:size(Y4,2)
    vs = Violin({Y4(:,i)},...
        i,...
        'HalfViolin','full',...% left, full
        'QuartileStyle','shadow',... % boxplot, none
        'DataStyle', 'scatter',... % scatter, none
        'ShowNotches', false,...
        'ShowMean', true,...
        'ShowMedian', true,...
        'MedianMarkerSize', 96,...
        'MarkerSize',80, ...
        'ViolinAlpha', {.3, .9},...
        'LineWidth', 3, ...
        'ViolinColor', {cc1(i, :)});
end
axis([0.5 2.5 0.5 1.3])
xticks([1 2])
% xlabel({'First line';'Second line'})
xticklabels({'Neighboring \newline illuminants','Opposite \newline illuminants'})
if relative
    ylabel('CCI relative to Control')
else
    ylabel('Color Constancy Index')

end


set(gca, 'FontSize', 34, 'fontname','Times New Roman');