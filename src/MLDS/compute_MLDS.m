function [targetCompetitorFit, logLikelyFit, predictedResponses] = ...-
    compute_MLDS(Data, competitorIndices, nTrials, colors, saveFig)

% Fit data for each target (column of Data) and then, plot the results. 
nTargets = size(Data, 2);
lengthLinY = max(competitorIndices(:));
targetCompetitorFit = cell(nTargets,1);
logLikelyFit = cell(nTargets,1);
predictedResponses = cell(nTargets,1);
for j = 1:nTargets
    nTrialsPerPair = nTrials(:,j); 
    theResponses = Data(:,j);   
    [targetCompetitorFit{j}, logLikelyFit{j}, predictedResponses{j}] = ...
     MLDSColorSelection(competitorIndices,theResponses,nTrialsPerPair, ...
     lengthLinY); 
    if saveFig
        MLDSColorSelectionPlot(competitorIndices,theResponses,...
        nTrialsPerPair,targetCompetitorFit{j},predictedResponses{j}, ...
        colors, saveFig);
    end
    clear theResponses nTrialsPerPair
end