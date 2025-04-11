% Load and prepare the data
% Assuming the data is loaded into a table called 'data'
% In practice, you would import your CSV file first with:
data = readtable('LocalSurround_Control_outdoor_ccip_colored.csv');

% First, let's examine the structure of our dataset
disp('Data Structure:');
disp('Number of samples:');
height(data)
disp('Variables in the dataset:');
data.Properties.VariableNames

% Convert categorical variables to factors
data.Participant = categorical(data.Participant);
data.LocalSurround = categorical(data.LocalSurround);
data.Baseline = categorical(data.Baseline);
data.Illuminant = categorical(data.Illuminant);
data.Direction = categorical(data.Direction);
data.Scene = categorical(data.Scene);

% Create a dataset containing only the rows with CCIdiff values (non-baseline conditions)
nonBaselineData = data(data.CCIdiff ~= 0, :);

% Get unique values for each factor
uniqueLocalSurround = unique(nonBaselineData.LocalSurround);
uniqueIlluminant = unique(nonBaselineData.Illuminant);
uniqueDirection = unique(nonBaselineData.Direction);


% LocalSurround effect
figure;
boxplot(nonBaselineData.CCIdiff, nonBaselineData.LocalSurround);
title('Effect of LocalSurround on CCIdiff');
xlabel('LocalSurround');
ylabel('CCIdiff');

% Illuminant effect
figure;
boxplot(nonBaselineData.CCIdiff, nonBaselineData.Illuminant);
title('Effect of Illuminant on CCIdiff');
xlabel('Illuminant');
ylabel('CCIdiff');

% Direction effect
figure;
boxplot(nonBaselineData.CCIdiff, nonBaselineData.Direction);
title('Effect of Direction on CCIdiff');
xlabel('Direction');
ylabel('CCIdiff');

% For LocalSurround x Illuminant interaction
figure;
uniqueLocalSurrounds = categories(nonBaselineData.LocalSurround);
uniqueIlluminants = categories(nonBaselineData.Illuminant);
hold on;

% Create markers and colors for different combinations
markers = {'o', 's', 'd', '^'};
colors = {'r', 'g', 'b', 'k', 'm', 'c'};
colors_local = [189,183,107;
    147,112,219;
    205,92,92;
    95,158,160;
    125,125,125;]./255;

colors_ill = [0.2316 0.2833 0.4579; %blue
    0.38 0.2441 0.1305;  %yellow
    0.233 0.3 0.2003;  %green
    0.4152 0.2114 0.3961].^.45;

% Plot mean CCIdiff for each LocalSurround x Illuminant combination
for i = 1:length(uniqueLocalSurrounds)
    for j = 1:length(uniqueIlluminants)
        % Get data for this combination
        idx = nonBaselineData.LocalSurround == uniqueLocalSurrounds{i} & ...
              nonBaselineData.Illuminant == uniqueIlluminants{j};
        
        if any(idx)
            meanVal = mean(nonBaselineData.CCIdiff(idx));
            plot(i, meanVal, 'o', 'MarkerFaceColor', colors_ill(j, :),...
                'MarkerSize', 50, ...
                 'DisplayName', [char(uniqueLocalSurrounds{i}), '-o', char(uniqueIlluminants{j})]);
            hold on;
        end
    end
end

% Add labels and legend
set(gca, 'XTick', 1:length(uniqueLocalSurrounds), 'XTickLabel', uniqueLocalSurrounds);
xlabel('Local Surround');
ylabel('Mean CCIdiff');
title('Interaction: LocalSurround x Illuminant');
legend('show', 'Location', 'best');
grid on;
hold off;

% Similarly for LocalSurround x Direction
figure;
uniqueDirections = categories(nonBaselineData.Direction);
hold on;

% Plot mean CCIdiff for each LocalSurround x Direction combination
for i = 1:length(uniqueLocalSurrounds)
    for j = 1:length(uniqueDirections)
        % Get data for this combination
        idx = nonBaselineData.LocalSurround == uniqueLocalSurrounds{i} & ...
              nonBaselineData.Direction == uniqueDirections{j};
        
        if any(idx)
            meanVal = mean(nonBaselineData.CCIdiff(idx));
            plot(i, meanVal, [colors{j}, markers{1}], 'MarkerSize', 10, ...
                 'DisplayName', [char(uniqueLocalSurrounds{i}), '-', char(uniqueDirections{j})]);
        end
    end
end

% Add labels and legend
set(gca, 'XTick', 1:length(uniqueLocalSurrounds), 'XTickLabel', uniqueLocalSurrounds);
xlabel('Local Surround');
ylabel('Mean CCIdiff');
title('Interaction: LocalSurround x Direction');
legend('show', 'Location', 'best');
grid on;
hold off;

% For Illuminant x Direction
figure;
hold on;

% Plot mean CCIdiff for each Illuminant x Direction combination
for i = 1:length(uniqueIlluminants)
    for j = 1:length(uniqueDirections)
        % Get data for this combination
        idx = nonBaselineData.Illuminant == uniqueIlluminants{i} & ...
              nonBaselineData.Direction == uniqueDirections{j};
        
        if any(idx)
            meanVal = mean(nonBaselineData.CCIdiff(idx));
            plot(i, meanVal, [colors{j}, markers{1}], 'MarkerSize', 10, ...
                 'DisplayName', [char(uniqueIlluminants{i}), '-', char(uniqueDirections{j})]);
        end
    end
end

% Add labels and legend
set(gca, 'XTick', 1:length(uniqueIlluminants), 'XTickLabel', uniqueIlluminants);
xlabel('Illuminant');
ylabel('Mean CCIdiff');
title('Interaction: Illuminant x Direction');
legend('show', 'Location', 'best');
grid on;
hold off;


figure;
gscatter(1:height(nonBaselineData), nonBaselineData.CCIdiff, ...
    [nonBaselineData.LocalSurround, nonBaselineData.Illuminant]);
title('Interaction: LocalSurround x Illuminant');
xlabel('Observation Index');
ylabel('CCIdiff');
legend('Location', 'Best');

% LocalSurround x Direction
figure;
gscatter(1:height(nonBaselineData), nonBaselineData.CCIdiff, ...
    [nonBaselineData.LocalSurround, nonBaselineData.Direction]);
title('Interaction: LocalSurround x Direction');
xlabel('Observation Index');
ylabel('CCIdiff');
legend('Location', 'Best');

% 8. Check model assumptions
figure;
plotResiduals(fullModel);

% 9. Run post-hoc tests for significant factors
% For example, if LocalSurround is significant:
disp('Post-hoc comparisons for LocalSurround:');
% Create all pairwise combinations of LocalSurround levels
surroundLevels = categories(nonBaselineData.LocalSurround);
numLevels = length(surroundLevels);
for i = 1:(numLevels-1)
    for j = (i+1):numLevels
        % Create a contrast vector
        contrastVec = zeros(1, numLevels);
        contrastVec(i) = 1;
        contrastVec(j) = -1;
        
        % Test the contrast
        contrastName = [char(surroundLevels(i)), ' vs ', char(surroundLevels(j))];
        try

            % Get data for each level
            data1 = nonBaselineData.CCIdiff...
                (nonBaselineData.LocalSurround == surroundLevels{i});
            data2 = nonBaselineData.CCIdiff...
                (nonBaselineData.LocalSurround == surroundLevels{j});

            % Perform t-test
            [h, p, ci, stats] = ttest2(data1, data2);
            fprintf('%s vs %s: t = %.4f, p = %.4f\n', ...
                char(surroundLevels{i}), char(surroundLevels{j}), stats.tstat, p);

            disp([contrastName, ': p = ', num2str(contrastResult.pValue)]);
        catch
            disp(['Error testing contrast for ', contrastName]);
        end
    end
end

% 10. Alternative model considering random slopes
% This allows the effect of predictors to vary by participant
randomSlopeModel = fitlme(nonBaselineData, 'CCIdiff ~ LocalSurround*Illuminant*Direction + (Direction|Participant)');
disp('Random slope model:');
disp(randomSlopeModel);

% Compare with the original random intercept model
disp('Comparing random intercept vs. random slope models:');
compareRandomModels = compare(fullModel, randomSlopeModel);
disp(compareRandomModels);

% Determine the importance of the intercept by evaluating the model fit criteria
disp('Intercept importance analysis:');
disp(['Full model with intercept AIC: ', num2str(fullModel.ModelCriterion.AIC)]);
disp(['Model without intercept AIC: ', num2str(noInterceptModel.ModelCriterion.AIC)]);
disp(['Full model with intercept BIC: ', num2str(fullModel.ModelCriterion.BIC)]);
disp(['Model without intercept BIC: ', num2str(noInterceptModel.ModelCriterion.BIC)]);

if fullModel.ModelCriterion.AIC < noInterceptModel.ModelCriterion.AIC
    disp('The model WITH intercept has lower AIC, suggesting the intercept is important');
else
    disp('The model WITHOUT intercept has lower AIC, suggesting the intercept might not be necessary');
end

% Save results
% In practice, you might want to save these results
% writetable(anovaResults, 'anova_results.csv');