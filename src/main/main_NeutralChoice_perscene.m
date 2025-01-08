clear
% close all
outdoor = 0;
relative = 1;
SZ = 1000;
competitors_lab =...
    [59.4233    2.7620   18.5559
    59.4768    1.3495    9.1942
    59.5303   -0.0630   -0.1675
    60.6783    0.8103  -13.5027
    61.8264    1.6836  -26.8379];

if outdoor==0

    tests_c = ...
        convertCharsToStrings({'control_khaki_indoor_ab.mat', ...
        'control_rose_indoor_ab.mat', ...
        'control_purple_indoor_ab.mat','control_teal_indoor_ab.mat'});

    tests = ...
        convertCharsToStrings({'khaki_indoor_ab.mat', ...
        'rose_indoor_ab.mat', ...
        'purple_indoor_ab.mat','teal_indoor_ab.mat'});

else

    tests_c = ...
        convertCharsToStrings({'control_khaki_outdoor_ab.mat', ...
        'control_rose_outdoor_ab.mat', ...
        'control_purple_outdoor_ab.mat','control_teal_outdoor_ab.mat'});

    tests = ...
        convertCharsToStrings({'khaki_outdoor_ab.mat', ...
        'rose_outdoor_ab.mat', ...
        'purple_outdoor_ab.mat','teal_outdoor_ab.mat'});

end



colors_ill_clockwise = {'Yellow', 'Red','Blue', 'Green'};
colors_ill_anticlockwise = {'Green', 'Yellow', 'Red','Blue'};
colors_ill_oppclockwise = {'Blue', 'Green','Yellow', 'Red'};
colors_ill_oppanticlockwise = {'Red', 'Blue', 'Green','Yellow'};

order_illuminants = {'Neutral', 'Blue', 'Yellow','Green', 'Red'};


local = {"Control Khaki", "Khaki", "Control Rose", "Rose",...
        "Control Purple", "Purple", "Control Teal","Teal"};


cci_control = cell(length(tests_c), 1);
cci_surround = cell(length(tests), 1);
for i=1:length(tests)
    cci_surround{i} = load(tests{i});
end
for i=1:length(tests_c)
    cci_control{i} = load(tests_c{i});
end


experiments_ill = [0.2316 0.2833 0.4579; %blue
    0.38 0.2441 0.1305;  %yellow
    0.233 0.3 0.2003;  %green
    0.4152 0.2114 0.3961; %red
    0.2912 0.2705 0.2614]; %neutral

%% Plot the local surround CCI and the control %%%%%%%%%%%%%%%%%%%%%%%%%%%%
condition_names = local;

c = [189,183,107;
    205,92,92;
    147,112,219;
    95,158,160;
    125,125,125]./255;

for i=1:length(cci_surround)
    aux = cci_surround{i}.new_lab_color(:, 1);
    data2{i} = cell2mat(aux);
end
for i=1:length(cci_control)
    aux = cci_control{i}.new_lab_color(:, 1);
    data1{i} = cell2mat(aux);
end


figure
%% Plot the b channel of the MATCH and 
%% the extrem competitors from the achromatic
plot([-1 15], competitors_lab(3,3).*ones(2,1), ...
    'Color',(2.*experiments_ill(end, :)).^.45,'LineStyle','--', ...
    'LineWidth', 6);hold on
plot([-1 15], competitors_lab(1,3).*ones(2,1), ...
    'Color',(2.*experiments_ill(2, :)).^.45,'LineStyle','--', ...
    'LineWidth', 6);hold on
plot([-1 15], competitors_lab(end,3).*ones(2,1), ...
    'Color',(2.*experiments_ill(1, :)).^.45,'LineStyle','--',...
    'LineWidth', 6);hold on

for i=1:length(cci_control)
    scatter((2*(i-1) + 2*(i-1)).*ones(length(data1{i}), 1), ...
       data1{i}(:, 3),SZ, ...
       'MarkerFaceColor',c(i, :), 'MarkerEdgeColor',[.15 .15 .15]);hold on
end

for i=1:length(cci_surround)
    scatter((2*(i-1) + 1 + 2*(i-1)).*ones(length(data2{i}), 1), ...
       data2{i}(:, 3),SZ, ...
       'MarkerFaceColor',c(i, :), 'MarkerEdgeColor',[.5 .5 .5]);hold on
end

for i=1:length(cci_control)
    for j=1:size(data1{i}, 1)
    plot([(2*(i-1) + 2*(i-1))...
          2*(i-1) + 1 + 2*(i-1)], ...
          [data1{i}(j, 3) data2{i}(j, 3)],'k--');hold on
    end
end

ylabel('b*');
aa = 1:length(cci_control);
num_xticks = sort([2*(aa-1)+2*(aa-1) 2*(aa-1)+2*(aa-1)+1]); 
xticks(num_xticks)
xticklabels(condition_names)
axis([2*(aa(1)-1)+2*(aa(1)-1)-1 2*(aa(end)-1)+2*(aa(end)-1)+2 ...
    -30 20])
set(gca, 'FontSize', 30, 'FontName','L M Roman10');


