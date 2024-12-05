function MLDSColorSelectionPlot(thePairs,theResponses,nTrialsPerPair,targetCompetitorFit,predictedResponses, colors, saveFig)
% function MLDSColorSelectionPlot(thePairs,theResponses,nTrialsPerPair,targetCompetitorFit,predictedResponses)
%
% Plots the inferred position of the target and all the competitors.
% If the theoretical probabilities are passed, they are plotted against the actual probabilities
% computed from the data.
%
%
% 5/3/12  dhb  Optional scatterplot of theory against measurements.
% 6/12/13 ar   Changed function names and added comments.
% 4/07/14 ar   Option to save the figures.
% 7/27/14 dhb  Update savefig to use more portable FigureSave.
global root_name participant illuminant

% Plot theoretical vs. actual probabilities?
if (nargin < 5 || isempty(predictedResponses))
    SCATTERPLOT = 0;
else
    SCATTERPLOT = 1;
end

% Display data and the fits.
theData = [thePairs theResponses./nTrialsPerPair];
fprintf('\n Target data set\n');
%disp(theData);
fprintf('\n Target fits\n');
%disp(targetCompetitorFit);

% Plot the inferred position of the target and the competitors.
f = figure; clf;
if (SCATTERPLOT)
    hold on
end

% axis([0 6 0 max(targetCompetitorFit)]);
for i = 2:length(targetCompetitorFit) % size(colors,1)
    %     plot(1:length(targetCompetitorFit(2:end)),targetCompetitorFit(2:end),'o','MarkerSize',10,'MarkerEdgeColor',[.2 .2 .2],'MarkerFaceColor',[.2 .2 .2].^(1/2));
    plot(i-1,targetCompetitorFit(i),'o','MarkerSize',10,'MarkerEdgeColor',colors(i-1, :),'MarkerFaceColor',colors(i-1, :).^(1/2.2)); hold on
end
plot(1:length(targetCompetitorFit(2:end)),targetCompetitorFit(1)*ones(size(targetCompetitorFit(2:end))),'k--','LineWidth',1.5);
title(sprintf('Inferred Target Position'));
xlabel('Competitor #');
ylabel('Representation');
axis('square');
set(gca,  'FontSize', 11, 'fontname','LM Mono Light 10');
saveas(gcf,strcat('mlds_', root_name, num2str(participant), '_ill', num2str(illuminant), '.svg'))
% export_fig(strcat('mlds_', root_name, num2str(participant), '_ill', num2str(illuminant)), '-eps')
close;

% Compute the probabilities from the data.
% Plot them vs. predicted probabilites from the fit.
if (SCATTERPLOT)
    figure;hold on
    theDataProb = theResponses./ nTrialsPerPair;
    plot(theDataProb,predictedResponses,'ro','MarkerSize',10); % ,'MarkerFaceColor','r'
    plot([0 1],[0 1],'k');
    axis('square')
    axis([0 1 0 1]);
    title(sprintf('Quality of MLDS model'));
    xlabel('Measured probabilities');
    ylabel('Predicted probabilities');
    set(gca,  'FontSize', 11, 'fontname', 'LM Mono Light 10');
    saveas(gcf,strcat('probabilities_', root_name, num2str(participant), '_ill', num2str(illuminant), '.svg'))
    %     export_fig(strcat('probabilities_', root_name, num2str(participant), '_ill', num2str(illuminant)), '-eps')
    close;
    
    figure;
    plot([log(targetCompetitorFit(3)-1e-4) log(targetCompetitorFit(end))], [2 2], 'Color', [150 150 150]./255);hold on
    for i = 3:length(targetCompetitorFit) % size(colors,1)
        %     plot(1:length(targetCompetitorFit(2:end)),targetCompetitorFit(2:end),'o','MarkerSize',10,'MarkerEdgeColor',[.2 .2 .2],'MarkerFaceColor',[.2 .2 .2].^(1/2));
        plot(log(targetCompetitorFit(i)), 2, 'o','MarkerSize',11,'MarkerEdgeColor',colors(i-1, :),'MarkerFaceColor',colors(i-1, :).^(1/2.2)); hold on
    end
    plot(log(targetCompetitorFit(1)), 2, 'o','MarkerSize',11,'MarkerEdgeColor',colors(i-1, :),'MarkerFaceColor','g');
    axis([log(targetCompetitorFit(3)-1e-4) log(targetCompetitorFit(end)) 1.9 2.1]); axis equal; axis off
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gca,  'FontSize', 11, 'fontname', 'LM Mono Light 10');
    saveas(gcf,strcat('logperceptual_', root_name, num2str(participant), '_ill', num2str(illuminant), '.svg'))
    %     export_fig(strcat('perceptual_', root_name, num2str(participant), '_ill', num2str(illuminant)), '-eps')
    close;
end
if saveFig
    FigureSave('MLDSOutput', f, 'pdf');
end
end