function new_lab_color = ...
    perceptual2Lab(targetCompetitorFit,targetCompetitorFit_match, ...
                   Colors_lab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if targetCompetitorFit_match > targetCompetitorFit(end)
    targetCompetitorFit_match = targetCompetitorFit(end);
end
new_lab_color = interp1(targetCompetitorFit,Colors_lab,...
targetCompetitorFit_match);



