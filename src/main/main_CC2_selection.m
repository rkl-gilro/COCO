%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%  Read and store the data from the csv files:
%%
%% Experiment|Trial|Red|Green|Blue|Xpos|Ypos|Zpos|Xtheta|Ytheta|Ztheta
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
folder0 = 'Data/CC2/';                                       % folder where the competitors are stored
folder = 'Data/CC2/selection_cubes_neutral_yellow/';
subdir = dir(folder);

calibration_file = ['/home/raquel/Documents/repositories/'...
    'CC2-main/rsrc/results_ColorCharacterization/'...
    'colorcharacteriation_mat/'...
    'Calib_UnrealStandard521_Varjo_13_06_2024_4_LUT_dE.mat'];

num_sessions =   1.*ones(7,1);
N = 5;                                                                     % num of competitors
competitorIndices = nchoosek(1:N, 2);
lab_111 = 1;                                                               % convert to LAB using default Matlab -1- or the actual calibration -0-
cci_space = 'Lab';                                                         % compute CCI in 'Lab' coordinates or only chromaticity coordinates AB
saveFig = [];
show = 1;                                                                  % set to 1 for plotting figures
C = 100;                                                                   % size of the scatter
save_file_name = strcat('testing_cubes_', cci_space,'.mat');

%% Neutral, blue, Yellow, Green and Red
experiments_pos = {'02_Neutral-Yellow-Control','02_Neutral-Yellow'};       % cubes
% experiments_pos = {'02_Neutral-Yellow-Spheres', ...                        % spheres
%     '02_Neutral-Yellow-Spheres-Control'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the illuminants
% Neutral-1-, Yellow-2-
experiments_ill = [0.2470    0.2163    0.2502;
    0.3508    0.2076    0.1267];
ill_neutral = [0.2470    0.2163    0.2502];

%% Give number labels to the experiments, same as illuminants
experiments_vals = [1, 2];

%% Read the files with the color sample competitors (reflectances)
%% In this case for both neutral and yellow illuminant the competitors 
%% are the same
T = readtable( strcat(folder0,'Yellow.csv') );
Competitors{1}  = cell2mat( table2cell(T) ) ;
Competitors{2}  = Competitors{1};

%% Read csv files and store in CSVtable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSVtable = cell(length(num_sessions), 1);
Part_names = cell(length(num_sessions), 1);
T = [];
for i = 1:length(num_sessions)
    for j = sum( num_sessions(1:i-1) )+1:sum( num_sessions(1:i) )
        T = [T; readtable([folder subdir(j+2).name])];
        subdir(j+2).name
    end
    %% Important: discard wrong trials!
    % r_color = find(T(:, 7) == 0);
    % T(r_color,:) = [];

    CSVtable{i} = T;

    aux_name = split(subdir(j+2).name, '_');
    Part_names{i} = aux_name{1};
    T = [];
end


%% Get only the combinations and the selection from participants data
Chosen_lizard_RGB = cell(length(CSVtable), 1);
Chosen_lizard_num = cell(length(CSVtable), 1);
Shown_lizard_RGB  = cell(length(CSVtable), 1);
Shown_lizard_num  = cell(length(CSVtable), 1);

Illuminants = cell(length(CSVtable), 1);

N_location_comp = zeros(N, N);
Kula_location = zeros(N, N);
Kula_location_achr = zeros(N, N);

Kula_loc_achr_blue = zeros(1, 3);

for i=1:length(CSVtable)

    %% Define illuminants (integers) based on experiments name
    Illuminants{i} = nan(size(CSVtable{i}, 1),1);
    for j = 1:length(experiments_vals)
        pos = find(strcmp(experiments_pos{j},CSVtable{i}.Time));
        Illuminants{i}(pos) = j;
    end
    %% Define lizard colours (integers) based on competitors values
    [Chosen_lizard_RGB{i}, Chosen_lizard_num{i}, ...
        Shown_lizard_RGB{i}, Shown_lizard_num{i}] = ...
        readorganise_csvdata0( CSVtable{i}, Competitors, experiments_pos);


    for k=1:N
        N_location_comp(:, k) = ...
            [N_location_comp(1,k)+...
            sum( (Chosen_lizard_num{i} == k).* ...
            (Shown_lizard_num{i}.S1 == k)' );
            N_location_comp(2,k)+...
            sum( (Chosen_lizard_num{i} == k).* ...
            (Shown_lizard_num{i}.S2 == k)' );
            N_location_comp(3,k)+...
            sum( (Chosen_lizard_num{i} == k).* ...
            (Shown_lizard_num{i}.S3 == k)' );
            N_location_comp(4,k)+...
            sum( (Chosen_lizard_num{i} == k).* ...
            (Shown_lizard_num{i}.S4 == k)' );
            N_location_comp(5,k)+...
            sum( (Chosen_lizard_num{i} == k).* ...
            (Shown_lizard_num{i}.S5 == k)' )];

        Kula_location(:, k) = ...
          [Kula_location(1,k)+sum( (Shown_lizard_num{i}.S1 == k) );...
           Kula_location(2,k)+sum( (Shown_lizard_num{i}.S2 == k) );...
           Kula_location(3,k)+sum( (Shown_lizard_num{i}.S3 == k) );...
           Kula_location(4,k)+sum( (Shown_lizard_num{i}.S4 == k) );...
           Kula_location(5,k)+sum( (Shown_lizard_num{i}.S5 == k) )];
        
        %% Particular configuration: blue under yellow illuminant and achr
        %% under netural illuminant. Choosing the blue
        Kula_loc_achr_blue(1, :) = [
            Kula_loc_achr_blue(1, 1)+...
        sum((Shown_lizard_num{i}.S1 == 4) & (Shown_lizard_num{i}.S5 == 3)) ...
             Kula_loc_achr_blue(1, 2)+sum(...
          (Shown_lizard_num{i}.S1 == 4) & (Shown_lizard_num{i}.S5 == 3)&...
          (Chosen_lizard_num{i} == 4)') ...
           Kula_loc_achr_blue(1, 3)+sum(...
          (Shown_lizard_num{i}.S1 == 4) & (Shown_lizard_num{i}.S5 == 3)&...
          (Chosen_lizard_num{i} == 3)')];

    end

end

group_names = {'Total', 'Selected'};
condition_names = {'Pos1', 'Pos2', 'Pos3', 'Pos4','Pos5'};
for i=1:N
    datap = [Kula_location(1,i), N_location_comp(1,i);...
        Kula_location(2,i) N_location_comp(2,i);...
        Kula_location(3,i) N_location_comp(3,i);...
        Kula_location(4,i) N_location_comp(4,i);...
        Kula_location(5,i) N_location_comp(5,i)];
    figure;
    % b = bar(positions, datap);
    % b(2).FaceColor = (2.*Competitors{1}(i, :)).^.35;
    % b(1).FaceColor = (2.*Competitors{1}(i, :)).^.45;
    c_aux = ...
        [(2.*Competitors{1}(i, :)).^.35;(Competitors{1}(i, :)).^.45];
    h = dabarplot(datap','group',[1,2],'xtlabels', condition_names,...
    'colors',c_aux,'scattersize',20,'bartype','stacked','numbers',2, ...
    'round',0); 
    ylabel('Frequency');
    legend([h.br(1,:)],group_names);    % add the legend manually
    xl = xlim; xlim([xl(1), xl(2)+1]);  % make more space for the legend
    ylim([0 190]); 
    set(gca, 'FontSize', 20, 'fontname','Times New Roman');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For each participant and illuminant compute Radonjic method (original)
%% Size for all variables # participants x # illuminants
Matrix              = cell(length(CSVtable), size(experiments_ill, 1));
targetCompetitorFit = cell(length(CSVtable), size(experiments_ill, 1));
logLikelyFit        = cell(length(CSVtable), size(experiments_ill, 1));
predictedResponses  = cell(length(CSVtable), size(experiments_ill, 1));
new_lab_color       = cell(length(CSVtable), size(experiments_ill, 1));

for i = 1:length(CSVtable)
    disp '------------------------------------------------------------'
    disp(['           Participant: ', num2str(i), '           '])

    for j = 1:size(experiments_ill, 1)

        pos = find( Illuminants{i} == j );
        M = zeros(N, N);
        for k = 1:length(pos)
            aux_rows = Chosen_lizard_num{i}(pos(k));
            aux_cols = [Shown_lizard_num{i}.S1(pos(k)) ...
                Shown_lizard_num{i}.S2(pos(k)) ...
                Shown_lizard_num{i}.S3(pos(k)) ...
                Shown_lizard_num{i}.S4(pos(k)) ...
                Shown_lizard_num{i}.S5(pos(k))];

            M(aux_rows, aux_cols) = M(aux_rows, aux_cols) + 1;
        end

        Matrix{i, j} = M;

        %% Rearrange M matrix into vectors for Radonjic algorithm
        count=1;
        for i_M = 1:N
            for j_M = i_M+1:N

                first_chosen(count) = M(i_M, j_M);
                total_times(count) = M(i_M, j_M) + M(j_M, i_M);
                count = count + 1;
            end
        end

        %% Radonjic algorithm to compute/estimate the winner/match under
        %% each illuminant
        total_times(total_times == 0) = 1;
        [targetCompetitorFit{i, j}, logLikelyFit{i, j},...
            predictedResponses{i, j}] = ...
            compute_MLDS(first_chosen', competitorIndices, total_times', ...
            Competitors{j}, saveFig);

        %% Convert the estimated locations of competitors and match in
        %% perceptual domain to LAB color space
        if lab_111
            competitors_lab = rgb2lab(Competitors{j}, ...
                'WhitePoint',[1 1 1], 'ColorSpace', 'linear-rgb');
            wp = [1 1 1];
        else
            XYZ_comp = applyCalibration_rgb2xyz(calibration_file, ...
            [Competitors{j};1 1 1]);
            competitors_lab = xyz2lab(XYZ_comp(1:end-1, :),...
            'WhitePoint', XYZ_comp(end, :));

            wp = XYZ_comp(end, :);
        end
        
        
        new_lab_color{i, j} = ...
            perceptual2Lab(targetCompetitorFit{i, j}{1}(2:end),...
            targetCompetitorFit{i, j}{1}(1), competitors_lab);% 
        
        %% Compute RGB values 
        if lab_111
            new_rgb_color{i, j} = ...
                lab2rgb(new_lab_color{i, j}, ...
                'WhitePoint',[1 1 1], 'ColorSpace', 'linear-rgb');
        else
            aux = lab2xyz(new_lab_color{i, j}, ...
                'WhitePoint', XYZ_comp(end, :));
            new_rgb_color{i, j} = ...
                applyCalibration_xyz2rgb(calibration_file, aux);
        end
    end
end


%% Compute the different color constancy indexes (CCI) and Brunswik ratio
%% Neutral, Blue, Yellow, Green and Red
% CCI = zeros(length(CSVtable), size(experiments_ill, 1));
figure;
for i =1:length(CSVtable)
    
    % no CCI under neutral illuminant
    for j = 2:size(experiments_ill, 1)
        values = [new_rgb_color{i, 1};Competitors{1}(3, :);...
            new_rgb_color{i, j};ill_neutral(1, :);...
            experiments_ill(j, :);1 1 1];
        XYZ = applyCalibration_rgb2xyz(calibration_file, values);

        %% In the reflected space:
        %% The match under neutral illuminant
        %% new_rgb_color{i, 1} .* experiments_ill(1, :)
        if lab_111
            neutralMatch_reflected = ...
                rgb2lab(new_rgb_color{i, 1}.*ill_neutral(1, :), ...
                'WhitePoint',[1 1 1],...
                'ColorSpace', 'linear-rgb');
        else
            xyz_aux = XYZ(1, :) .* XYZ(4, :);
            xyz_aux = xyz_aux ./ XYZ(4, 2);                                % divide by luminance of the illuminant (https://en.wikipedia.org/wiki/CIE_1931_color_space)

            neutralMatch_reflected = ...
                xyz2lab(xyz_aux,'WhitePoint', XYZ(end, :));                % the achromatic match under the neutral illuminant

        end


        %% Competitors{1}(3/5, :) .* experiments_ill(1, :)                 % the achromatic under the neutral illuminant
        if lab_111
            neutralTest_reflected = ...
            rgb2lab(Competitors{1}(3, :).*ill_neutral(1, :), ...
            'WhitePoint',[1 1 1],...
            'ColorSpace', 'linear-rgb');
        else
            xyz_aux = XYZ(2, :) .* XYZ(4, :);
            xyz_aux = xyz_aux ./ XYZ(4, 2);

            neutralTest_reflected = ...
                xyz2lab(xyz_aux, 'WhitePoint', XYZ(end, :));
        end
        
        
        %% new_rgb_color{i, j} .* experiments_ill(j, :)                    % the match under the colored illuminant
        if lab_111
            coloredMatch_reflected = ...
                rgb2lab(new_rgb_color{i, j}.*experiments_ill(j, :), ...
                'WhitePoint',[1 1 1],...
                'ColorSpace', 'linear-rgb');
        else
            xyz_aux = XYZ(3, :) .* XYZ(5, :);
            xyz_aux = xyz_aux ./ XYZ(5, 2);

            coloredMatch_reflected = ...
                xyz2lab(xyz_aux, 'WhitePoint', XYZ(end, :));
        end
        
        
        %% Competitors{1}(5, :) .* experiments_ill(j, :)                   % the achromatic under the colored illuminant
        if lab_111
            coloredTest_reflected = ...
                rgb2lab(Competitors{1}(3, :).*experiments_ill(j, :), ...
                'WhitePoint',[1 1 1],...
                'ColorSpace', 'linear-rgb');
        else
            xyz_aux = XYZ(2, :) .* XYZ(5, :);
            xyz_aux = xyz_aux ./ XYZ(5, 2);

            coloredTest_reflected = ...
                xyz2lab(xyz_aux, 'WhitePoint', XYZ(end, :));
        end


        if new_lab_color{i, 1} ~= -1 & new_lab_color{i, j}  ~= -1
            
            if strcmp(cci_space, 'Lab')
                dist_a  = ...
                    norm(coloredTest_reflected-neutralTest_reflected);
                dist_a1 = ...
                    norm(coloredTest_reflected-neutralMatch_reflected);    % the neutral match instead of the test
                dist_b  = ...
                    norm(coloredTest_reflected-coloredMatch_reflected);
                dist_c  = ...
                    norm(neutralTest_reflected-coloredMatch_reflected);
                dist_c1 = ...
                    norm(neutralMatch_reflected-coloredMatch_reflected);   % the neutral match instead of the test
            else
                dist_a  = norm(coloredTest_reflected(2:end)...
                    -neutralTest_reflected(2:end));
                dist_a1 = norm(coloredTest_reflected(2:end)...
                    -neutralMatch_reflected(2:end));                       % the neutral match instead of the test
                dist_b  = norm(coloredTest_reflected(2:end)...
                    -coloredMatch_reflected(2:end));
                dist_c  = norm(neutralTest_reflected(2:end)...
                    -coloredMatch_reflected(2:end));
                dist_c1 = norm(neutralMatch_reflected(2:end)...
                    -coloredMatch_reflected(2:end));                       % the neutral match instead of the test
            end

            CCI.cci(i, j)         = 1 - dist_b/dist_a;
            CCI.cci_neutral(i, j) = 1 - dist_b/dist_a1;
            CCI.BR(i, j)          = dist_c / dist_a;
            CCI.BR_neutral(i, j)  = dist_c1 / dist_a1;

            scatter3(coloredTest_reflected(2), ...
                coloredTest_reflected(3), ...
                coloredTest_reflected(1), C, ...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[0 .75 .75], ...
                'Marker', 'square');hold on

            scatter3(neutralTest_reflected(2), ...
                neutralTest_reflected(3), ...
                neutralTest_reflected(1), C, ...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[0.75 .75 .75], ...
                'Marker', 'square');hold on

            scatter3(coloredMatch_reflected(2), ...
                coloredMatch_reflected(3), ...
                coloredMatch_reflected(1), C, ...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[0 .75 .75]);hold on

            scatter3(neutralMatch_reflected(2), ...
                neutralMatch_reflected(3), ...
                neutralMatch_reflected(1), C, ...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[0.75 .75 .75]);hold on

            xlabel('a');ylabel('b');zlabel('L');
            set(gca, 'FontSize', 20, 'fontname','Times New Roman');
        else
            'Error: CCI not computed!!'
            CCI.cci(i, j) = -1;
            CCI.cci_neutral(i, j) = -1;
            CCI.BR(i, j) = -1;
            CCI.BR_neutral(i, j) = -1;
        end

    end
end


%% Save in a mat file
save(save_file_name, 'Part_names', 'new_lab_color', 'CCI',...
    'targetCompetitorFit', 'logLikelyFit','predictedResponses', 'Matrix');
