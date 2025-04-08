% Mixed Model Analysis Pipeline for Experimental Data with Incomplete Design
% This script provides a comprehensive analysis of participant data where not all
% participants completed all conditions, using mixed models and ANOVA.

function data = analyzeExperimentData(csvFilePaths, dependentVars)
    % Main function to analyze experimental data from multiple CSV files
    %
    % Parameters:
    % csvFilePaths - Cell array of paths to CSV files containing participant data
    % dependentVars - Cell array of dependent variables to analyze (e.g., {'CCI', 'MatchPosition'})
    
    % Check input arguments
    if nargin < 2
        dependentVars = {'CCI', 'MatchPosition'};
        fprintf('Using default dependent variables: CCI and MatchPosition\n');
    end
    
    % Step 1: Data Loading and Preprocessing
    fprintf('==== Step 1: Loading and Preprocessing Data ====\n');
    data0 = loadAndPreprocessData(csvFilePaths);
    data = removevars(data0,{'ColorShiftControl'});

    % Step 2: Explore data structure
    fprintf('\n==== Step 2: Data Structure Exploration ====\n');
    exploreDataStructure(data);
    
    % Step 3: Run mixed models for each dependent variable
    fprintf('\n==== Step 3: Mixed Models Analysis ====\n');
    for i = 1:length(dependentVars)
        fprintf('\n----- Analyzing %s -----\n', dependentVars{i});
        runMixedModels(data, dependentVars{i});
    end
    
    % Step 4: Visualize key relationships
    fprintf('\n==== Step 4: Visualizing Key Relationships ====\n');
    for i = 1:length(dependentVars)
        visualizeResults(data, dependentVars{i});
    end
    
    % Step 5: Save processed data and results
    fprintf('\n==== Step 5: Saving Results ====\n');
    saveResults(data);
    
    fprintf('\nAnalysis complete!\n');
end

function combinedData = loadAndPreprocessData(csvFilePaths)
    % Loads and preprocesses data from multiple CSV files
    
    % Initialize empty table
    combinedData = table();
    
    % Load and combine all files
    fprintf('Loading %d CSV files...\n', length(csvFilePaths));
    for i = 1:length(csvFilePaths)
        [~, filename, ~] = fileparts(csvFilePaths{i});
        fprintf('Processing file %d/%d: %s\n', i, length(csvFilePaths), filename);
        
        try
            % Read the CSV file
            fileData = readtable(csvFilePaths{i});
            
            % Check and handle expected columns
            expectedCols = {'CCI', 'Participant', 'LocalSurround', 'Illuminant', 'Direction', ...
                           'ColorShift', 'ColorShiftLS', 'MatchPosition', ...
                           'Baseline', 'Scene'};
            
            missingCols = setdiff(expectedCols, fileData.Properties.VariableNames);
            if ~isempty(missingCols)
                warning('File %s is missing columns: %s', filename, strjoin(missingCols, ', '));
                % Add missing columns with NaN values
                for j = 1:length(missingCols)
                    fileData.(missingCols{j}) = NaN(height(fileData), 1);
                end
            end
            
            % Append file data to combined data
            combinedData = [combinedData; fileData];
            
        catch e
            warning('Error processing file %s: %s', filename, e.message);
        end
    end
    
    % Convert categorical variables to categorical type
    categoricalVars = {'Participant', 'LocalSurround', 'Illuminant', 'Direction', 'Baseline', 'Scene'};
    for i = 1:length(categoricalVars)
        if ismember(categoricalVars{i}, combinedData.Properties.VariableNames)
            % Check if the column contains non-numeric data or should be treated as categorical
            if iscell(combinedData.(categoricalVars{i})) || ...
               ischar(combinedData.(categoricalVars{i})) || ...
               isa(combinedData.(categoricalVars{i}), 'categorical') || ...
               length(unique(combinedData.(categoricalVars{i}))) < 10
                
                combinedData.(categoricalVars{i}) = categorical(combinedData.(categoricalVars{i}));
                fprintf('Converted %s to categorical with %d levels\n', ...
                    categoricalVars{i}, length(categories(combinedData.(categoricalVars{i}))));
            end
        end
    end
    
    % Check for missing values
    missingCount = sum(ismissing(combinedData));
    if any(missingCount > 0)
        fprintf('\nMissing values detected:\n');
        for i = 1:length(combinedData.Properties.VariableNames)
            if missingCount(i) > 0
                fprintf('  %s: %d missing values (%.1f%%)\n', ...
                    combinedData.Properties.VariableNames{i}, ...
                    missingCount(i), ...
                    (missingCount(i)/height(combinedData))*100);
            end
        end
    else
        fprintf('No missing values detected.\n');
    end
    
    % Remove rows with NaN in key variables
    keyVars = {'Participant', 'CCI', 'MatchPosition'};
    keyVarsPresent = intersect(keyVars, combinedData.Properties.VariableNames);
    nanRows = any(ismissing(combinedData(:, keyVarsPresent)), 2);
    if any(nanRows)
        fprintf('Removing %d rows with missing values in key variables\n', sum(nanRows));
        combinedData = combinedData(~nanRows, :);
    end
    
    fprintf('Final dataset: %d observations from %d participants\n', ...
        height(combinedData), length(unique(combinedData.Participant)));
end

function exploreDataStructure(data)
    % Explores the structure of the data, checking available conditions and design balance
    
    fprintf('Dataset dimensions: %d rows x %d columns\n', height(data), width(data));
    
    % Check design balance
    fprintf('\nChecking design balance:\n');
    
    % Look at key factors
    factors = {'Direction', 'LocalSurround', 'Illuminant', 'Baseline', 'Scene'};
    availableFactors = intersect(factors, data.Properties.VariableNames);
    
    for i = 1:length(availableFactors)
        factor = availableFactors{i};
        levels = categories(data.(factor));
        if isempty(levels)
            levels = unique(data.(factor));
        end
        
        fprintf('%s has %d levels: %s\n', factor, length(levels), strjoin(string(levels), ', '));
        
        % Check how many participants completed each level
        for j = 1:length(levels)
            participantsInLevel = unique(data.Participant(data.(factor) == levels{j}));
            fprintf('  Level %s: %d participants\n', string(levels{j}), length(participantsInLevel));
        end
    end
    
    % Check crossover between factors (for main factors only)
    mainFactors = intersect({'Direction', 'LocalSurround', 'Illuminant'}, availableFactors);
    if length(mainFactors) >= 2
        fprintf('\nCross-tabulation of main factors:\n');
        for i = 1:length(mainFactors)-1
            for j = i+1:length(mainFactors)
                crosstab_result = crosstab(data.(mainFactors{i}), data.(mainFactors{j}));
                fprintf('Crosstab of %s x %s:\n', mainFactors{i}, mainFactors{j});
                disp(crosstab_result);
            end
        end
    end
    
    % Check distribution of dependent variables
    dependentVars = {'CCI', 'MatchPosition'};
    availableDVs = intersect(dependentVars, data.Properties.VariableNames);
    
    fprintf('\nDependent variable statistics:\n');
    for i = 1:length(availableDVs)
        dv = availableDVs{i};
        fprintf('%s: Mean = %.3f, SD = %.3f, Range = [%.3f to %.3f]\n', ...
            dv, mean(data.(dv), 'omitnan'), std(data.(dv), 'omitnan'), ...
            min(data.(dv)), max(data.(dv)));
        
        % Test normality
        try
            [h, p] = lillietest(data.(dv));
            fprintf('  Lilliefors normality test: p = %.4f, %s\n', p, ...
                ternaryOp(h==0, 'Normal distribution', 'Non-normal distribution'));
        catch
            fprintf('  Could not perform normality test on %s\n', dv);
        end
    end
end

function runMixedModels(data, dependentVar)
    % Runs mixed models analysis for the specified dependent variable
    
    % Check if fitlme is available (requires Statistics and Machine Learning Toolbox)
    if ~exist('fitlme', 'file')
        error(['Mixed models require the Statistics and Machine Learning Toolbox. ', ...
               'Please install it or use alternative methods.']);
    end
    
    % Check if the dependent variable exists
    if ~ismember(dependentVar, data.Properties.VariableNames)
        error('Dependent variable %s not found in the dataset', dependentVar);
    end
    
    % Identify available fixed factors
    availableFactors = intersect({'Direction', 'LocalSurround', 'Illuminant', 'Baseline', 'Scene'}, ...
                                data.Properties.VariableNames);
    
    % Build formulas with increasing complexity
    fprintf('Fitting mixed models for %s with Participant as random effect...\n', dependentVar);
    
    % Model 1: Individual main effects
    fprintf('\n--- Individual main effects ---\n');
    modelResults = struct();
    
    for i = 1:length(availableFactors)
        factor = availableFactors{i};
        formula = sprintf('%s ~ 1 + %s + (1|Participant)', dependentVar, factor);
        fprintf('Model formula: %s\n', formula);
        
        try
            lme = fitlme(data, formula);
            modelResults.(factor) = lme;
            
            % Display results
            fprintf('AIC: %.2f, BIC: %.2f\n', lme.ModelCriterion.AIC, lme.ModelCriterion.BIC);
            fprintf('Fixed Effects:\n');
            disp(lme.Coefficients(2:end, [1,2,3,6]));  % Display estimate, SE, tStat, pValue
            
            fprintf('\nANOVA table (Type III tests):\n');
            anovalme = anova(lme, 'DFMethod', 'Satterthwaite');
            disp(anovalme);

            % Test significance of the main effect using anova
            reduced = fitlme(data, sprintf('%s ~ 1 + (1|Participant)', dependentVar));
            [~, pValue, ~, ~] = compare(reduced, lme, 'CheckNesting', true);
            fprintf('Main effect of %s: p = %.4f\n', factor, pValue);

        catch e
            fprintf('Error fitting model for %s: %s\n', factor, e.message);
        end
    end
    
     
    % Model 3: Full model with all main effects (and selected interactions if significant)
    fprintf('\n--- Full model with all main effects ---\n');
    
    % Build formula with all main effects
    fullFormula = sprintf('%s ~ 1', dependentVar);
    for i = 1:length(availableFactors)
        fullFormula = [fullFormula, ' + ', availableFactors{i}];
    end
    fullFormula = [fullFormula, ' + (1|Participant)'];
    
    try
        fullModel = fitlme(data, fullFormula);
        fprintf('Full model formula: %s\n', fullFormula);
        fprintf('AIC: %.2f, BIC: %.2f\n', fullModel.ModelCriterion.AIC, fullModel.ModelCriterion.BIC);
        fprintf('Fixed Effects:\n');
        disp(fullModel.Coefficients(2:end, [1,2,3,6]));
        
        % Perform Type III tests for significance of each factor
        fprintf('\nANOVA table (Type III tests):\n');
        anovaResults = anova(fullModel, 'DFMethod', 'Satterthwaite');
        disp(anovaResults);
        
    catch e
        fprintf('Error fitting full model: %s\n', e.message);
        
        % Try a more conservative model if full model fails
        fprintf('Trying simplified model...\n');
        try
            % Use only the first 3 factors or fewer if less are available
            numFactors = min(3, length(availableFactors));
            simplifiedFormula = sprintf('%s ~ 1', dependentVar);
            for i = 1:numFactors
                simplifiedFormula = [simplifiedFormula, ' + ', availableFactors{i}];
            end
            simplifiedFormula = [simplifiedFormula, ' + (1|Participant)'];
            
            simplifiedModel = fitlme(data, simplifiedFormula);
            fprintf('Simplified model formula: %s\n', simplifiedFormula);
            fprintf('Fixed Effects:\n');
            disp(simplifiedModel.Coefficients(2:end, [1,2,3,6]));
            
            % Perform Type III tests
            fprintf('\nANOVA table (Type III tests):\n');
            anovaResults = anova(simplifiedModel, 'DFMethod', 'Residual'); % , 'DFMethod', 'Satterthwaite'
            disp(anovaResults);
            
        catch e2
            fprintf('Error fitting simplified model: %s\n', e2.message);
        end
    end
end

function visualizeResults(data, dependentVar)
    % Creates visualizations for key relationships
    colors = [189,183,107;
    147,112,219;
    205,92,92;
    95,158,160;
    125,125,125;]./255;

    experiments_ill = [0.2316 0.2833 0.4579; %blue
    0.38 0.2441 0.1305;  %yellow
    0.233 0.3 0.2003;  %green
    0.4152 0.2114 0.3961].^.45;

    % Check if the dependent variable exists
    if ~ismember(dependentVar, data.Properties.VariableNames)
        fprintf('Cannot visualize: %s not found in dataset\n', dependentVar);
        return;
    end
    
    % Identify available fixed factors
    availableFactors = intersect({'Direction', 'LocalSurround', 'Illuminant', 'Baseline', 'Scene'}, ...
                                data.Properties.VariableNames);
    colorFactors{1} = [.3 .3 .3;.6 .6 .6];% baseline
    colorFactors{2} = [.2 .4 .6;.3 .5 .8;.2 .4 .6]; % direction
    colorFactors{3} = experiments_ill;%illuminant
    colorFactors{4} = colors;%local surround
    colorFactors{5} = [.6 .5 .6;.2 .9 .6];%scene

    if isempty(availableFactors)
        fprintf('No factors available for visualization\n');
        return;
    end
    
    fprintf('Creating visualizations for %s...\n', dependentVar);
    
    % 1. Main effects plots
    for i = 1:length(availableFactors)
        factor = availableFactors{i};
        plotMainEffect(data, dependentVar, factor, colorFactors{i});
    end
    
    % 2. Interaction plots for primary factors (if available)
    primaryFactors = intersect({'Direction', 'LocalSurround', 'Illuminant'}, availableFactors);
    colorFactors2{1} = [.2 .4 .6;.3 .5 .8]; % direction
    colorFactors2{2} = experiments_ill;%illuminant
    colorFactors2{3} = colors;%local surround
    if length(primaryFactors) >= 2
        for i = 1:length(primaryFactors)-1
            for j = i+1:length(primaryFactors)
                plotInteraction(data, dependentVar, primaryFactors{i}, ...
                    primaryFactors{j}, colorFactors2{j});
            end
        end
    end
    
    % 3. Faceted plots for main factor combinations
    if length(primaryFactors) >= 2
        plotFacetedResults(data, dependentVar, primaryFactors, colors);
    end
end

function plotMainEffect(data, dependentVar, factor, colors)
    % Plots the main effect of a factor on the dependent variable
    
    % Create a new figure
    figure('Name', sprintf('Main Effect of %s on %s', factor, dependentVar));
    
    % Calculate means and standard errors
    levels = categories(data.(factor));
    if isempty(levels)
        levels = unique(data.(factor));
    end
    
    means = zeros(length(levels), 1);
    sems = zeros(length(levels), 1);
    ns = zeros(length(levels), 1);
    
    for i = 1:length(levels)
        levelData = data.(dependentVar)(data.(factor) == levels{i});
        levelData_cell{i} = levelData;
        means(i) = mean(levelData, 'omitnan');
        sems(i) = std(levelData, 'omitnan') / sqrt(sum(~isnan(levelData)));
        ns(i) = sum(~isnan(levelData));
    end
    
    % Create bar plot with error bars
    h = dabarplot(levelData_cell,'xtlabels', levels,...
    'colors',colors,'errorbars','SE', 'scattersize', 40);

    % bar(means,'FaceColor',colors,'EdgeColor',[.3 .3 .3]);
    % hold on;
    % errorbar(1:length(levels), means, sems,'k', 'LineStyle', 'none');
    
    % Add labels and formatting
    title(sprintf('Effect of %s on %s', factor, dependentVar));
    xlabel(factor);
    ylabel(dependentVar);
    xticks(1:length(levels));
    xticklabels(levels);
    xtickangle(45);
    
    % Add sample size labels
    for i = 1:length(levels)
        text(i, means(i) + sems(i) + 0.05 * range(means), ...
             sprintf('n=%d', ns(i)), ...
             'HorizontalAlignment', 'center');
    end
    
    grid on;
    box off;
    hold off;
    set(gca, 'FontSize', 40, 'fontname','L M Roman10');
    % Save figure
    saveas(gcf, sprintf('main_effect_%s_on_%s.png', factor, dependentVar));
end

function plotInteraction(data, dependentVar, factor1, factor2, colors)
    % Plots the interaction between two factors on the dependent variable
    
    % Create a new figure
    figure('Name', sprintf('Interaction of %s and %s on %s', factor1, factor2, dependentVar));
    
    % Get levels of each factor
    levels1 = categories(data.(factor1));
    if isempty(levels1)
        levels1 = unique(data.(factor1));
    end
    
    levels2 = categories(data.(factor2));
    if isempty(levels2)
        levels2 = unique(data.(factor2));
    end
    
    % Calculate means and standard errors
    means = zeros(length(levels1), length(levels2));
    sems = zeros(length(levels1), length(levels2));
    
    for i = 1:length(levels1)
        for j = 1:length(levels2)
            cellData = data.(dependentVar)(data.(factor1) == levels1{i} & data.(factor2) == levels2{j});
            means(i,j) = mean(cellData, 'omitnan');
            sems(i,j) = std(cellData, 'omitnan') / sqrt(sum(~isnan(cellData)));
        end
    end
    
    % Create interaction plot
    for j = 1:length(levels2)
        errorbar(1:length(levels1), means(:,j), sems(:,j), '-s', 'Markersize', 10,'LineWidth',...
            1.5,'MarkerEdgeColor',colors(j, :), 'MarkerFaceColor',colors(j, :), 'Color',colors(j, :));
        hold on;
    end
    
    % Add labels and formatting
    title(sprintf('Interaction of %s and %s on %s', factor1, factor2, dependentVar));
    xlabel(factor1);
    ylabel(dependentVar);
    xticks(1:length(levels1));
    xticklabels(levels1);
    xtickangle(45);
    legend(cellstr(levels2), 'Location', 'best');
    grid on;
    box off;
    hold off;
    
    % Save figure
    saveas(gcf, sprintf('interaction_%s_%s_on_%s.png', factor1, factor2, dependentVar));
end

function plotFacetedResults(data, dependentVar, factors, colors)
    % Creates a faceted plot showing relationships between multiple factors
    
    % Limit to 3 factors maximum for visualization
    factors = factors(1:min(3, length(factors)));
    
    % Create a new figure
    figure('Name', sprintf('Faceted Plot of %s', dependentVar), 'Position', [100, 100, 1000, 600]);
    
    % Get levels of each factor
    levels = cell(length(factors), 1);
    for i = 1:length(factors)
        levels{i} = categories(data.(factors{i}));
        if isempty(levels{i})
            levels{i} = unique(data.(factors{i}));
        end
    end
    
    % Determine subplot layout based on number of factors
    switch length(factors)
        case 1
            rows = 1; cols = 1;
        case 2
            rows = length(levels{2}); cols = 1;
        case 3
            rows = length(levels{2}); cols = length(levels{3});
    end
    
    % Create subplots
    spIdx = 1;
    
    for i = 1:rows
        for j = 1:cols
            subplot(rows, cols, spIdx);
            
            % Filter data based on factor levels for this subplot
            if length(factors) == 1
                subplotData = data;
                titleStr = dependentVar;
            elseif length(factors) == 2
                subplotData = data(data.(factors{2}) == levels{2}{i}, :);
                titleStr = sprintf('%s = %s', factors{2}, string(levels{2}{i}));
            else
                subplotData = data(data.(factors{2}) == levels{2}{i} & ...
                                  data.(factors{3}) == levels{3}{j}, :);
                titleStr = sprintf('%s = %s, %s = %s', ...
                                 factors{2}, string(levels{2}{i}), ...
                                 factors{3}, string(levels{3}{j}));
            end
            
            % Create boxplot of factor 1
            boxplot(subplotData.(dependentVar), subplotData.(factors{1}));
            
            % Add title and labels
            title(titleStr);
            xlabel(factors{1});
            ylabel(dependentVar);
            
            % Rotate x labels if needed
            if length(levels{1}) > 3
                xtickangle(45);
            end
            
            spIdx = spIdx + 1;
        end
    end
    
    % Add overall title
    sgtitle(sprintf('Effects on %s by %s', dependentVar, strjoin(factors, ', ')));
    
    % Adjust spacing
    tight_layout = true;
    if tight_layout && exist('tight_subplot', 'file')
        % Use tight_subplot if available
        tight_subplot(rows, cols, 'Padding', 0.05, 'MarginLeft', 0.1, 'MarginRight', 0.05);
    else
        % Otherwise use standard spacing adjustment
        pos = get(gcf, 'Position');
        set(gcf, 'Position', [pos(1), pos(2), pos(3), pos(4)]);
    end
    
    % Save figure
    saveas(gcf, sprintf('faceted_plot_%s.png', dependentVar));
end

function saveResults(data)
    % Saves the processed data and results
    
    % Save the processed data
    writetable(data, 'processed_data.csv');
    fprintf('Processed data saved to processed_data.csv\n');

end

% Helper function for ternary operation
function result = ternaryOp(condition, trueVal, falseVal)
    if condition
        result = trueVal;
    else
        result = falseVal;
    end
end