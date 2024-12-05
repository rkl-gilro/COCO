function CI = Color_Index4( ill1, ill2, ill_r1, ill_r2 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% CCI
%% Color Constancy Index (CCI; Ling & Hurlbert, 2008)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compute the CI metric
Sc = ill_r2 - ill2;
Sm = ill_r1 - ill1;
Sp = ill2 - ill1;

d = 1;
thr = 10^(-3);
if norm(Sc) < norm(Sm) && dot(Sc, Sm)<thr
    d = -1
end

% CI =  1 - (dot(Sc-Sm, -Sp) / norm( Sp ));
CI =  1 - (norm(Sc-Sm)*d / norm( Sp ));