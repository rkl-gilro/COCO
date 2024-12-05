function CI = Color_Index3( ill1, ill2, ill_r1, ill_r2 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% BR constancy index
%% Brunswick ratio (BR; Hansen et al., 2007; Smithson & Zaidi, 2004; Yang & Shevell, 2002)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compute the CI metric
CI =  norm( ill_r2 - ill_r1 ) / norm( ill2 - ill1 );
