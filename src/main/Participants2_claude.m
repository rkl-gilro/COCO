% Comprehensive Correlation Analysis for Participants

% Load data
data = readtable('LocalSurround_Control_outdoor_BRpn_colored.csv');

% Convert categorical variables
data.Participant = categorical(data.Participant);
data.Illuminant = categorical(data.Illuminant);
data.Baseline = categorical(data.Baseline);

% Get unique participants and illuminants
participants = unique(data.Participant);
illuminants = unique(data.Illuminant);

% Preallocate correlation matrices
intraParticipantCorr = nan(length(participants), length(illuminants));
interParticipantCorr = nan(length(participants), length(participants));

% Intra-Participant Correlation
% Correlation of color shifts across different illuminants for each participant
for p = 1:length(participants)
    participant_data = data(data.Participant == participants(p), :);
    
    % Prepare color shifts matrix
    color_shifts = nan(length(illuminants), 1);
    
    % Calculate mean color shift for each illuminant
    for i = 1:length(illuminants)
        illuminant_data = participant_data(participant_data.Illuminant == illuminants(i), :);
        
        % Check if data exists for this illuminant
        if ~isempty(illuminant_data)
            color_shifts(i) = nanmean(illuminant_data.CCIdiff);
        end
    end
    
    % Calculate correlation if enough data
    if sum(~isnan(color_shifts)) > 1
        % Remove NaN values before correlation calculation
        valid_shifts = color_shifts(~isnan(color_shifts));
        valid_corr = corrcoef(valid_shifts);
        
        % Store mean correlation
        intraParticipantCorr(p, 1:size(valid_corr,1)) = valid_corr;
    end
end

% Inter-Participant Correlation
% Correlation of participants' color shift patterns
for p1 = 1:length(participants)
    for p2 = 1:length(participants)
        % Get aggregated data for both participants
        p1_data = data(data.Participant == participants(p1) & ...
            data.Baseline == 'Suppression', :);
        p2_data = data(data.Participant == participants(p2) & ...
            data.Baseline == 'Suppression', :);
        
        % Aggregate color shifts by illuminant
        p1_agg = grpstats(p1_data, 'LocalSurround', {'mean'}, 'DataVars', 'CCIdiff');
        p2_agg = grpstats(p2_data, 'LocalSurround', {'mean'}, 'DataVars', 'CCIdiff');
        % p1_agg2 = grpstats(p1_data, 'Illuminant', {'mean'}, 'DataVars', 'CCIdiff');
        % p2_agg2 = grpstats(p2_data, 'Illuminant', {'mean'}, 'DataVars', 'CCIdiff');
        % Ensure both have same illuminants
        [C,ia,ib] = intersect(p1_agg.LocalSurround, p2_agg.LocalSurround);
        % [C2,ia2,ib2] = intersect(p1_agg2.Illuminant, p2_agg2.Illuminant);
        
        % Calculate correlation if enough data
        if length(C) > 1
            interParticipantCorr(p1, p2) = ...
                corr(p1_data.CCI,...
            p2_data.CCI);
            %     corr(p1_agg.mean_CCIdiff(ia),...
            % p2_agg.mean_CCIdiff(ib));
            % corr([p1_agg.mean_CCIdiff(ia);p1_agg2.mean_CCIdiff(ia2)],...
            % [p2_agg.mean_CCIdiff(ib);p2_agg2.mean_CCIdiff(ib2)]);
        end
    end
end

% Clean up NaN values
intraParticipantCorr(isnan(intraParticipantCorr)) = 0;
interParticipantCorr(isnan(interParticipantCorr)) = 0;

% Visualization of Correlation Matrices
% Intra-Participant Correlation Heatmap
% figure;
% heatmap(cellstr(string(participants)), cellstr(string(illuminants)), intraParticipantCorr);
% title('Intra-Participant Correlation Across Illuminants');
% xlabel('Participants');
% ylabel('Illuminants');

% Inter-Participant Correlation Heatmap
figure;
heatmap(cellstr(string(participants)), cellstr(string(participants)), interParticipantCorr);
title('Inter-Participant Correlation');
xlabel('Participants');
ylabel('Participants');
colormap("parula");

% Summary Statistics
disp('Intra-Participant Correlation Summary:');
disp(mean(intraParticipantCorr, 'all'));
disp('Standard Deviation of Intra-Participant Correlation:');
disp(std(intraParticipantCorr, 0, 'all'));

disp('Inter-Participant Correlation Summary:');
disp(mean(interParticipantCorr, 'all'));
disp('Standard Deviation of Inter-Participant Correlation:');
disp(std(interParticipantCorr, 0, 'all'));

% Export correlation matrices
writematrix(intraParticipantCorr, 'intra_participant_correlation.csv');
writematrix(interParticipantCorr, 'inter_participant_correlation.csv');