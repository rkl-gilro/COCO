%% response1 = MLDSSimulateResponse(x,y1,y2,sigma,mapFunction)
%
% Simulates a trial given a target and a competitor pair.
% The passed mapFunction simulates the effect of context change
% between x domain and y domain
function response = MLDSSimulateResponse(x,y1,y2,sigma,mapFunction)

yOfX = mapFunction(x);
diff1 = y1-yOfX;
diff2 = y2-yOfX;
if (abs(diff1)-abs(diff2) + normrnd(0,sigma) <= 0)
    response = 1;
else
    response = 0;
end

end





