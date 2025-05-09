
function [targetCompetitorFit, logLikelyFit, predictedResponses] = MLDSColorSelection_sigma(thePairs,theResponses,nTrialsPerPair, numberOfCompetitors)
% function [targetCompetitorFit, logLikelyFit, predictedResponses] = MLDSColorSelection(thePairs,theResponses,nTrialsPerPair, numberOfCompetitors)
%
% Main MLDS function. Takes as input set of observer's responses and returns the inferred position of the target and each competitors.
%
% Input:
%   thePairs -            competitor pairs.
%   theResponses -        set of responses for this pair (number of times first
%                         competitor is chosen).
%   nTrialsPerPair -      total number of trials run. Vector of same size as theResponses.
%   numberOfCompetitors - number of competitors in the set.
%
% Output:
%   targetCompetitorFit - inferred positions for target and competitors.
%                         the target is the first entry and then come the competitors.
%   logLikelyFit -        log likelihood of the fit.
%   predictedResponses -  responses predicted from the fit.
%
% NOTE.  A number of hard coded parameters in this routine were developed for
% the specific color selection experiments we have been doing over the past
% year (2012-2013).  It may be that this routine is brittle with respect to
% major changes in the experimental design, and fussing with these parameters
% might help:
%  sigma
%  sigmaFactor
%  trySpacings
%
% 5/3/12  dhb  Capture predicted probabilities and pass them along.
% 6/13/13 ar   Added comments.
% 2/12/15 dhb, ar Enforce sigma spacing for first found competitor position.
% 2/23/15 dhb, ar More comments.  And add an input arg error check.

%% The number of responses in theResponses cannot ever exceed nTrialsPerPair.
% Check this and throw error message if it does not hold.
if (length(theResponses(:)) ~= length(nTrialsPerPair(:)))
    error('Passed theResponses and nTrialsPerPair must be of same length');
end
if (any(theResponses > nTrialsPerPair))
    error('An entry of input theResponses exceeds passed nTrialsPerPair.  No good!');
end

%% Set fixed parameters
%
% Standard deviation for the MLDS solution.
% This determines the scale of the solution.
sigma = 0.05; % .1; % 0.05;
disp(['Sigma: ', num2str(sigma)])
% Determine minimum size of interval between
% solution elements, relative to sigma.  That
% is, the minimum spacing will be sigma/sigmaFactor.
sigmaFactor = 4;

%% Set up parameters for search.

% Enforce constraint that competitor portion of
% the solution is monotonic with
% respect to the nominal target positions, with
% a reasonable separation (in sigma units).
%
% For the second through last competitor, we
% do this using the linear constraint feature
% that fmincon is set up to use, by generating
% a matrix A that takes the differences between
% adjacent entries of the solution vector.
% By requiring that Ax < b, we enforce that
% each difference in the solution is at
% least -b.  That is, the first entry of
% A computes C2-C3, which then enforces
% C2-C3 < -b => C3-C2 > b, which is what
% we want.
%
% Because the target is the first entry of
% the solution, it has no constraint.  We
% build the constraint only for the differences
% between the competitor solutions.
% A = zeros(numberOfCompetitors-2,numberOfCompetitors+1);
% for i = 2:numberOfCompetitors-2
%     A(i,i+1) = 1;
%     A(i,i+2) = -1;
% end
%
% % This is the minimum interval size for use with the
% % the A matrix above.
% b =  zeros(numberOfCompetitors-2,1); % -sigma/sigmaFactor*ones(numberOfCompetitors-2,1); %  -.02 .* ones(numberOfCompetitors-2,1); %

% Set spacings for initializing the search.
% Try different ones in the hope that we thus avoid local minima in
% the search.
%
% Note that these spacings are hard coded and were determined
% from experience.  Thus they are rather
% specific to the type of color selection experiment we
% have been doing.  It is possible that there would be a
% cleverer thing to do here.
% trySpacings = [0.5 1 2];

% Standard fmincon options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','sqp');

disp 'log instead of log10!!!!!! (MLDSComputeLogLikelihood)'
%% Search
%
% We search over various initial spacings and take the best result.
% There are two loops.  One sets the positions of the competitors
% in the solution.  The other tries different places for the target
% in the solution.
%
% In terms of variable names below, Y indicates competitor positoins,
% while X is (unfortunately
% maxLogLikely = -Inf;
% for k1 = 1:length(trySpacings)
%     % Choose initial competitor positions based on current spacing to try.
%     initialCompetitorPositions = linspace(0,trySpacings(k1)*numberOfCompetitors*sigma,numberOfCompetitors);
%
%     % Choose a set of initial target positions for the current initial competitor positions.
%     tryTargetPositions(k1,:) = linspace(min(initialCompetitorPositions),max(initialCompetitorPositions),2*numberOfCompetitors);
%
%     % We don't actually search on the position of the first competitor, since the model predictions don't change
%     % if we add a constant to all of the positions.  Thus we pull out the first competitor from the parameters
%     % vector.  Saving it here just allows us to plop it back in for computations later.
%     y1 = initialCompetitorPositions(1);
%
%     % Loop over all initial target positions that we are going to try
%     for k = 1:length(tryTargetPositions(k1,:))

% Initialize the parameters vector for the search, by prepending the initial target position
% to the list of all but the first initial competitor positions.
 % tryTargetPositions(k1,k);
L = 1;
initialCompetitorPositions = 0:L/(numberOfCompetitors-1):L;
initialTargetPosition = initialCompetitorPositions( round(length(initialCompetitorPositions)/2) ); % round(length(initialCompetitorPositions)/2)

y1 = initialCompetitorPositions(1);
yn = initialCompetitorPositions(end);

%% Target position, the competitors position (except first and final), and the sigma

initialParams = [initialCompetitorPositions(1) initialCompetitorPositions(2:end-1) 0.1];
% initialParams = [initialCompetitorPositions(1) initialCompetitorPositions(2:end) 0.1];

A = zeros( length(initialParams)-3, length(initialParams));
for i = 1: length(initialParams)-3
    A(i,i+1) =  1;
    A(i,i+2) = -1;
end

% This is the minimum interval size for use with the
% the A matrix above.
b =  -(1e-1).*ones( length(initialParams)-3,1);

%         initialTargetPosition = 3;
%         initialParams = [initialTargetPosition 1:5 0.1];
%
% Get reasonable upper and lower bound. These are most easily computed from the initial parameters.
% We enforce that the competitor solutions head off in the positive direction, but the target can
% be anywhere.  (We take 100 times the maximum value in the intial parameters to equal 'anyware');
%
% Because the first competitor is at 0, the rest cannot be lower
% than sigma/sigmaFactor.  This has the effect of enforcing C2 >
% sigma/sigmaFactor.  And since C3 > C2 etc, using this as a lower
% bound for C3 etc is OK.
epsilon = (1e-1);
vlb = zeros(size(initialParams))    + epsilon; % (sigma/sigmaFactor)*ones(size(initialParams));
vub = L.*ones(size(initialParams))  - epsilon; % 100*max(abs(initialParams))*ones(size(initialParams)); %100.*ones(size(initialParams)); %

% Remember that the searched vector is [T C2 C3 ...].  We don't
% want a strong lower bound constraint on the target, so we make
% this the negative our our large upper bound.  This lets the
% target, in effect, go anywhere.  That's what we want.
%         vlb(1) = -vub(1);


maxLogLikely = -Inf; 
% for i=1:length(initialCompetitorPositions)
    

initialParams = [initialTargetPosition initialCompetitorPositions(2:end-1) 0.1];
% initialParams = [initialTargetPosition initialCompetitorPositions(2:end) 0.1];
    
% Run the search
fitParams = fmincon(@(x)FitContextFCScalingFun_sigma(x,y1,thePairs,theResponses,nTrialsPerPair, yn),initialParams,A,b,[],[],vlb,vub,[],options);

% Extract target and competetior positions from the solution, and prepend the first competitor postion to
% final positions (see comment above).
fitTargetPosition = fitParams(1);
fitCompetitorPositions = [y1 fitParams(2:end-1) yn];
% fitCompetitorPositions = [y1 fitParams(2:end-1)];
sigma = fitParams(end);

% Compute log likelihood for this solution.  We need this so that we can keep track of the best
% solution that comes out of the multiple starting points.
% Save this solution if it's better than the current best.
%         temp = MLDSComputeLogLikelihood_sigma(thePairs,theResponses,nTrialsPerPair,fitTargetPosition,fitCompetitorPositions,sigma);

[logLikelyFit,predictedResponses] = MLDSComputeLogLikelihood_sigma(thePairs,theResponses,nTrialsPerPair,fitTargetPosition,fitCompetitorPositions,sigma);
targetCompetitorFit = [fitTargetPosition, fitCompetitorPositions];
disp(['Sigma: ', num2str(sigma)])

%         if (logLikelyFit > maxLogLikely)
%             maxLogLikely = logLikelyFit;
% %             [logLikelyFit,predictedResponses] = MLDSComputeLogLikelihood_sigma(thePairs,theResponses,nTrialsPerPair,fitTargetPosition,fitCompetitorPositions,sigma);
%             targetCompetitorFit = [fitTargetPosition, fitCompetitorPositions];
%         end
% %     end
% end

end

function f = FitContextFCScalingFun_sigma(x,y1,thePairs,theResponses,nTrials, L)
%function f = FitContextFCScalingFun(x,y1,thePairs,theResponses,nTrials,sigma)
% The error function we are minimizing in the numerical search.
% Computes the negative log likelyhood of the current solution i.e. inferred
% positio of the target and the competitors.
% Input:
%   x           - current target position fit.
%   y1          - current competitor position fits.
%   thePairs    - competitor pairs.
%   theResponses- set of responses for this pair (number of times first
%                 competitor is chosen.
%   nTrialsPerPair - total number of trials run.
%   sigma          - fixed standard deviation
% Output:
%   f - negative log likelihood for the current solution.


% Sanity check.
if (any(isnan(x)))
    error('Entry of x is NaN');
end

% Reorganize the solution vector.
% Pull the target from the solution, put the fixed competitor 1 back in the competitor array.
xFit = x(1);
yFit = [y1 x(2:end-1) L];
sigma = x(end);

A = zeros( length(yFit)-1, length(yFit));
for i = 1: length(yFit)-1
    A(i,i) =  -1;
    A(i,i+1) = 1;
end

ordered = A * yFit';
ordered(ordered > 0) = 1; 
if sigma <= 0 || ( ~isempty( ordered(ordered < 0) ) && sum(ordered(ordered < 0)) < 0 )
    f = Inf;
    return 
end
% compute negative log likelyhood of the current solution
logLikely = MLDSComputeLogLikelihood_sigma(thePairs,theResponses,nTrials,xFit,yFit,sigma);

f = -logLikely;

end
