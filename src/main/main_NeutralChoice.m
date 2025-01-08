clear
% close all
outdoor = 0;
control = 0;
relative = 1;
SZ = 600;
competitors_lab =...
   [59.4233    2.7620   18.5559
   59.4768    1.3495    9.1942
   59.5303   -0.0630   -0.1675
   60.6783    0.8103  -13.5027
   61.8264    1.6836  -26.8379];

if outdoor==0

    if control 
        tests = ...
            convertCharsToStrings({'control_khaki_indoor_ab.mat', ...
            'control_rose_indoor_ab.mat', ...
            'control_purple_indoor_ab.mat','control_teal_indoor_ab.mat'});
    else
        tests = ...
            convertCharsToStrings({'khaki_indoor_ab.mat', ...
            'rose_indoor_ab.mat', ...
            'purple_indoor_ab.mat','teal_indoor_ab.mat'});
    end
else
    if control
    tests = ...
        convertCharsToStrings({'control_khaki_outdoor_ab.mat', ...
        'control_rose_outdoor_ab.mat', ...
        'control_purple_outdoor_ab.mat','control_teal_outdoor_ab.mat'});
    else
        tests = ...
        convertCharsToStrings({'khaki_outdoor_ab.mat', ...
        'rose_outdoor_ab.mat', ...
        'purple_outdoor_ab.mat','teal_outdoor_ab.mat'});
    end
end



colors_ill_clockwise = {'Yellow', 'Red','Blue', 'Green'};
colors_ill_anticlockwise = {'Green', 'Yellow', 'Red','Blue'};
colors_ill_oppclockwise = {'Blue', 'Green','Yellow', 'Red'};
colors_ill_oppanticlockwise = {'Red', 'Blue', 'Green','Yellow'};

order_illuminants = {'Neutral', 'Blue', 'Yellow','Green', 'Red'};

if control
    local = {"Control Khaki", "Control Rose", ...
        "Control Purple", "Control Teal"};
else
    local = {"Khaki", "Rose", ...
        "Purple", "Teal"};
end

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
    0.4152 0.2114 0.3961; %red
   0.2912 0.2705 0.2614]; %neutral

%% Plot the local surround CCI and the control %%%%%%%%%%%%%%%%%%%%%%%%%%%%
condition_names = local;
condition_names2 = {"Opposing", "Neighboring"};

c = [189,183,107;
    205,92,92;
    147,112,219;
    95,158,160;
    125,125,125]./255;

for i=1:length(cci_surround)
    % data1{i} = mean(cci_surround{i}.CCI.cci_neutral(:, 2:end)', 2);
    aux = cci_surround{i}.new_lab_color(:, 1);
    data1{i} = cell2mat(aux); % sum( abs([59.5303   -0.0630   -0.1675] - cell2mat(aux)), 2 )./3;
end

yellow_lab = sum( abs([59.5303   -0.0630   -0.1675] - ...
    [59.4233    2.7620   18.5559]), 2 )./3;      % yellow
blue_lab = sum( abs([59.5303   -0.0630   -0.1675] - ...
    [61.8264    1.6836  -26.8379]), 2 )./3;      % blue

figure
plot([0 5], competitors_lab(3,3).*ones(2,1), ...
    'Color',experiments_ill(end, :).^.45,'LineStyle','--', ...
    'LineWidth', 4);hold on
plot([0 5], competitors_lab(1,3).*ones(2,1), ...
    'Color',experiments_ill(2, :).^.45,'LineStyle','--', ...
    'LineWidth', 4);hold on
plot([0 5], competitors_lab(end,3).*ones(2,1), ...
    'Color',experiments_ill(1, :).^.45,'LineStyle','--',...
    'LineWidth', 4);hold on
for i=1:length(cci_surround)
    scatter(i.*ones(length(data1{i}), 1), data1{i}(:, 3),SZ, ...
        'MarkerFaceColor',c(i, :), 'MarkerEdgeColor',[.5 .5 .5]);hold on
end
ylabel('b*');
xticks([ 1 2 3 4])
xticklabels(condition_names)
set(gca, 'FontSize', 30, 'FontName','L M Roman10');


