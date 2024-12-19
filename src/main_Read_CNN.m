%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Read and store the data from the csv files:
%%  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
folder0 = 'Data/LocalSurround/';                                           % folder where the competitors are stored
folder = 'Data/Achromatic_Control_20241212/csv/';
subdir = dir(folder);
allFileNames = {subdir.name};

illuminants = {'Neutral', 'Blue', 'Yellow', 'Green', 'Red'};


N = 5;                                                                     % num of competitors
                                                 % convert to LAB using default Matlab -1- or the actual calibration -0-                                                     % compute CCI in 'Lab' coordinates or only chromaticity coordinates AB
C = 200; 
reflectance_space = 0;
save_file_name = strcat('control_indoor_achromatic.mat');
saveFig = [];
show = 1;                                                                  % set to 1 for plotting figures

%% Define the illuminants
% Neutral-1-, Blue-2-, Yellow-3-, Green-4- and Red-5-
experiments_ill = [0.2912    0.2705    0.2614;
    0.2316    0.2833    0.4579;
    0.3800    0.2441    0.1305;
    0.2330    0.3000    0.2003;
    0.4152    0.2114    0.3961];


CSVtable = cell(length(illuminants), 1);
T = [];
numFilesProcessed = 0;	% For fun, let's keep track of how many files we processed.
for i = 1:length(illuminants)
    pattern = illuminants(i);
    for k = 1 : length(allFileNames)
    	% Get this filename.
    	thisFileName = fullfile(subdir(k).folder, allFileNames{k});
    	% See if it contains our required pattern.
    	if ~contains(thisFileName, pattern, 'IgnoreCase', true)
    		% Skip this file because the filename does not contain the required pattern.
    		continue;
        end
        T = [T; readtable(thisFileName)];
    	
    end
    CSVtable{i} = T;
    Lab_ill{i} = xyz2lab([T.x T.y T.z],'WhitePoint',[1 1 1]);
    aux_val = T.Value;
    aux_val(aux_val == 0) = 1;
    % RGB_ill{i} = hsv2rgb([T.Hue T.Chroma./aux_val T.Value./aux_val]);

    % H_prime = T.Hue ./60;
    % X = T.Chroma .* (1 - abs(mod(H_prime, 2) -1));
    % 
    % if( (0 <= H_prime) & (H_prime < 1))
    %     RGB1 = [T.Chroma X zeros(length(T.Chroma), 1)];
    % elseif(1 <= H_prime & H_prime < 2)
    %     RGB1 = [X T.Chroma zeros(length(T.Chroma), 1)];
    % elseif(2 <= H_prime & H_prime < 3)
    %     RGB1 = [zeros(length(T.Chroma), 1) T.Chroma X];
    % elseif(3 <= H_prime & H_prime < 4)
    %     RGB1 = [zeros(length(T.Chroma), 1) X T.Chroma];
    % elseif(4 <= H_prime & H_prime < 5)
    %     RGB1 = [X zeros(length(T.Chroma), 1) T.Chroma];
    % elseif(5 <= H_prime & H_prime < 6)
    %     RGB1 = [T.Chroma zeros(length(T.Chroma), 1) X];
    % end
    % 
    % m = T.Value - T.Chroma;
    % RGB = RGB1 + m;

    T = [];
end

achromatic = ...
rgb2lab([.145, .142, .145], 'Colorspace', 'linear-rgb', ...
    'WhitePoint', [1 1 1]);

achromatic2 = ...
rgb2lab([0.146317	0.145001	0.146794], 'Colorspace', 'linear-rgb', ...
    'WhitePoint', [1 1 1]);

if show
    figure
    for i = 1:length(illuminants)
        s = scatter3(CSVtable{i}, 'a', 'b', 'l');hold on
        s.SizeData = C;
        s.MarkerFaceColor = experiments_ill(i, :).^.45;
        s.MarkerEdgeColor = experiments_ill(i, :);
        s.MarkerFaceAlpha = .7;
        % , C, 'MarkerFaceColor',experiments_ill(i, :) 
        mean_val(i, :) = mean([CSVtable{i}.a CSVtable{i}.b CSVtable{i}.l]);
        ss = scatter3(mean_val(i, 1), mean_val(i, 2), ...
            mean_val(i, 3));hold on
        ss.SizeData = C*3;
        ss.MarkerFaceColor = experiments_ill(i, :).^.45;
        ss.MarkerEdgeColor = experiments_ill(i, :);

    end
    scatter3(achromatic(2), achromatic(3), achromatic(1), ...
        C*2,'filled','black');

    figure
    for i = 1:length(illuminants)
        s = scatter3(Lab_ill{i}(:, 2), Lab_ill{i}(:, 3), ...
            Lab_ill{i}(:, 1));hold on
        s.SizeData = C;
        s.MarkerFaceColor = experiments_ill(i, :).^.45;
        s.MarkerEdgeColor = experiments_ill(i, :);
        s.MarkerFaceAlpha = .7;

        mean_val2(i, :) = mean([Lab_ill{i}(:, 2), Lab_ill{i}(:, 3), ...
            Lab_ill{i}(:, 1)]);
        ss = scatter3(mean_val2(i, 1), mean_val2(i, 2), ...
            mean_val2(i, 3));hold on
        ss.SizeData = C*3;
        ss.MarkerFaceColor = experiments_ill(i, :).^.45;
        ss.MarkerEdgeColor = experiments_ill(i, :);

    end
    scatter3(achromatic(2), achromatic(3), achromatic(1), ...
        C*2,'filled','black');

end

aux = repmat(achromatic,5, 1);
deltaE00(aux', [aux(:, 1) mean_val(:, 1) mean_val(:, 2)]')
deltaE00(aux', [aux(:, 1) mean_val2(:, 1) mean_val2(:, 2)]')

save(save_file_name, 'achromatic', 'mean_val2', 'mean_val','CSVtable');
