%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Read and store the data from the csv files:
%%  Experiment|Trial|Red|Green|Blue|Xpos|Ypos|Zpos|Xtheta|Ytheta|Ztheta
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
folder0 = 'Data/LocalSurround/';                                           % folder where the competitors are stored
folder = 'Data/LocalSurround/Indoor/Teal/';
subdir = dir(folder);

calibration_file = ['/home/raquel/Documents/repositories/CC2-main/'...
    'rsrc/results_ColorCharacterization/'...
    'Calibration_UnrealStandard_HTCVive_14_03_2023_LUT_dE.mat'];

num_sessions =   [2 2 1 2 1]; % 2.*ones(7,1);%
N = 5;                                                                     % num of competitors
competitorIndices = nchoosek(1:N, 2);
lab_111 = 1;                                                               % convert to LAB using default Matlab -1- or the actual calibration -0-
cci_space = 'ab';                                                          % compute CCI in 'Lab' coordinates or only chromaticity coordinates AB
C = 200; 
reflectance_space = 0;
save_file_name = strcat('teal_indoor_',cci_space,'.mat');
saveFig = [];
show = 1;                                                                  % set to 1 for plotting figures


%% Define the illuminants
% Neutral-1-, Blue-2-, Yellow-3-, Green-4- and Red-5-
experiments_ill = [0.2912    0.2705    0.2614;
    0.2316    0.2833    0.4579;
    0.3800    0.2441    0.1305;
    0.2330    0.3000    0.2003;
    0.4152    0.2114    0.3961];

%% Neutral, blue, Yellow, Green and Red
experiments_pos = {'0.2912 0.2705 0.2614';
    '0.2316 0.2833 0.4579';
    '0.38 0.2441 0.1305';
    '0.233 0.3 0.2003';
    '0.4152 0.2114 0.3961'};

%% Give number labels to the experiments, same as illuminants
experiments_vals = [1, 2, 3, 4, 5];

%% Read the files with the color sample competitors (reflectances)
T = readtable( strcat(folder0,'Blue.csv') );
Competitors{2} = cell2mat( table2cell(T) ) ;
T = readtable( strcat(folder0,'Green.csv') );
Competitors{4}  = cell2mat( table2cell(T) ) ;
T = readtable( strcat(folder0,'Yellow.csv') );
Competitors{3}  = cell2mat( table2cell(T) ) ;
T = readtable( strcat(folder0,'Red.csv') );
Competitors{5}  =  cell2mat( table2cell(T) ) ;
T = readtable( strcat(folder0,'Neutral.csv') );
Competitors{1}  = cell2mat( table2cell(T) );

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
    r_color = find(T.Red == 0);
    r_lizard = find(strcmp(T.LizardName, ''));
    T([r_color; r_lizard],:) = [];

    CSVtable{i} = T;

    aux_name = split(subdir(j+2).name, '-');
    Part_names{i} = aux_name{1};

    T = [];
end


%% Get only the combinations and the selection from participants data
Chosen_lizard_RGB = cell(length(CSVtable), 1);
Chosen_lizard_num = cell(length(CSVtable), 1);
Shown_lizard_RGB  = cell(length(CSVtable), 1);
Shown_lizard_num  = cell(length(CSVtable), 1);

Illuminants = cell(length(CSVtable), 1);
for i=1:length(CSVtable)

    %% Define illuminants (integers) based on experiments name
    Illuminants{i} = nan(size(CSVtable{i}, 1),1);
    for j = 1:length(experiments_vals)
        pos = find(strcmp(experiments_pos{j},CSVtable{i}.Experiment));
        Illuminants{i}(pos) = j;
    end
    %% Define lizard colours (integers) based on competitors values
    [Chosen_lizard_RGB{i}, Chosen_lizard_num{i}, ...
      Shown_lizard_RGB{i}, Shown_lizard_num{i}] = ...
        readorganise_csvdata( CSVtable{i}, Competitors, experiments_pos);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For each participant and illuminant compute Radonjic method (original)
%% Size for all variables # participants x # illuminants
Matrix              = cell(length(CSVtable), size(experiments_ill, 1));
targetCompetitorFit = cell(length(CSVtable), size(experiments_ill, 1));
logLikelyFit        = cell(length(CSVtable), size(experiments_ill, 1));
predictedResponses  = cell(length(CSVtable), size(experiments_ill, 1));
new_lab_color       = cell(length(CSVtable), size(experiments_ill, 1));
Rsq1                = cell(length(CSVtable), size(experiments_ill, 1));

for i = 1:length(CSVtable)
    disp '------------------------------------------------------------'
        disp(['           Participant: ', num2str(i), '           '])
    for j = 1:size(experiments_ill, 1)

        pos = find( Illuminants{i} == j );
        
        M = zeros(N, N);                                                   % matrix to store the number of times the row was chosen against the column

        for k = 1:length(pos)
            aux_rows = Chosen_lizard_num{i}(pos(k));
            aux_cols = [Shown_lizard_num{i}.S1(pos(k)) ...
                Shown_lizard_num{i}.S2(pos(k)) ... 
                Shown_lizard_num{i}.S3(pos(k)) ...
                Shown_lizard_num{i}.S4(pos(k)) ...
                Shown_lizard_num{i}.S5(pos(k))];

            M(aux_rows, aux_cols) = M(aux_rows, aux_cols) + 1;

        end
        
        Matrix{i}{j} = M;
        
        %% Rearrange M matrix into vectors for Radonjic algorithm
        count=1;
        for i_M = 1:N
            for j_M = i_M+1:N
                
                first_chosen(count) = M(i_M, j_M);
                total_times(count) = M(i_M, j_M) + M(j_M, i_M);
                count = count + 1;
            end
        end

        disp(['---------- illuminant: ', num2str(j)])
        %% Radonjic algorithm to compute/estimate the winner/match under
        %% each illuminant
        total_times(total_times == 0) = 1;
        [targetCompetitorFit{i, j}, logLikelyFit{i, j},...
            predictedResponses{i, j}] = ...
    compute_MLDS(first_chosen', competitorIndices, total_times', ...
    Competitors{j}, saveFig);
        
        %% Compute model confidence, Pearson coeff
        Rsq = corrcoef( [(first_chosen ./ total_times)' ...
            cell2mat(predictedResponses{i, j})'] );

        if size(Rsq,1) < 2
            Rsq1{i, j} = Rsq;
        else
            Rsq1{i, j} = Rsq(2, 1);
        end
        
        if show
            figure(3)
            subplot(length(CSVtable), size(experiments_ill, 1), ...
                (i-1)*size(experiments_ill, 1) + j)
            plot(cell2mat(predictedResponses{i, j}), ...
                first_chosen ./ total_times, 'ro', 'MarkerSize',10);hold on
            xlabel('Data');ylabel('Predicted data');
            plot([0 1],[0 1],'k');
            axis('square')
            axis([0 1 0 1]);
            title(['Ill:', num2str(j),', R^2=', num2str(Rsq1{i, j})]);
        end

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
    targetCompetitorFit{i, j}{1}(1), competitors_lab);
        
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
if show
    figure;
end
for i =1:length(CSVtable)
    
    for j = 2:size(experiments_ill, 1)                                     % no CCI under neutral illuminant
        
        values = [new_rgb_color{i, 1};Competitors{1}(3, :);...
            new_rgb_color{i, j};experiments_ill(1, :);...
            experiments_ill(j, :);1 1 1];
        XYZ = applyCalibration_rgb2xyz(calibration_file, values);

        %% In the reflected space:
        %% The match under neutral illuminant
        %% new_rgb_color{i, 1} .* experiments_ill(1, :)
        if lab_111
            neutralMatch_reflected = ...
                rgb2lab(new_rgb_color{i, 1}.*experiments_ill(1, :), ...
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
            rgb2lab(Competitors{1}(3, :).*experiments_ill(1, :), ...
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
            
            if show
                scatter3(coloredTest_reflected(2), ...
                    coloredTest_reflected(3), ...
                    coloredTest_reflected(1), C, ...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',experiments_ill(j, :).^.45, ...
                    'Marker', 'square');hold on

                scatter3(neutralTest_reflected(2), ...
                    neutralTest_reflected(3), ...
                    neutralTest_reflected(1), C, ...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[.75 .75 .75], ...
                    'Marker', 'square');hold on

                scatter3(coloredMatch_reflected(2), ...
                    coloredMatch_reflected(3), ...
                    coloredMatch_reflected(1), C, ...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',experiments_ill(j, :).^.45);hold on

                scatter3(neutralMatch_reflected(2), ...
                    neutralMatch_reflected(3), ...
                    neutralMatch_reflected(1), C, ...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[0.75 .75 .75]);hold on

                xlabel('a');ylabel('b');zlabel('L');
                set(gca, 'FontSize', 20, 'fontname','Times New Roman');
            end
        else
            'Error: CCI not computed!!'
            CCI.cci(i, j) = -1;
            CCI.cci_neutral(i, j) = -1;
            CCI.BR(i, j) = -1;
            CCI.BR_neutral(i, j) = -1;
        end
        
    end
end

%% Save in a mat file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(save_file_name, 'Part_names', 'new_lab_color', 'CCI', 'Rsq1', ...
    'targetCompetitorFit','logLikelyFit','predictedResponses', 'Matrix');
