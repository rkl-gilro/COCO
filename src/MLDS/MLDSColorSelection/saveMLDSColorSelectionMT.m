function saveMLDSColorSelectionMT(competitorIndices,responses,numbertrials, nstimuli,fl)

  [targetCompetitorFit, logLikelyFit, predictedResponses] = MLDSColorSelection(competitorIndices,responses,numbertrials, nstimuli);
        
         save(fl)