function CI = Color_Index5_3D( ill1, ill2, ill_r1, ill_r2 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


AB = ill2 - ill1;
AP1 = ill_r1 - ill1;
AP2 = ill_r2 - ill1;
ill_r1_p = ill1 + dot(AP1,AB)/dot(AB, AB) * AB;
ill_r2_p = ill1 + dot(AP2,AB)/dot(AB, AB) * AB;

%% Compute the CI metric
disp('Penalize bigger than 1')
if norm( ill_r1_p - ill_r2_p ) > norm( ill1 - ill2 )
    CI = 2 - norm( ill_r1_p - ill_r2_p ) / norm( ill1 - ill2 );
else
    CI = norm( ill_r1_p - ill_r2_p ) / norm( ill1 - ill2 );
end
% CI = norm( ill_r1_p - ill_r2_p ) / norm( ill1 - ill2 );
    