function CI = Color_Index2( ill1, ill2, ill_r )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compute the CI metric
CI = 1 - ( norm( ill2 - ill_r ) / norm( ill2 - ill1 ) );


% figure;
% plot([ill1(1) ill2(1)], [ill1(2) ill2(2)], 'co-');hold on
% plot(ill_r(1), ill_r(2), '*k');hold on
% plot(x, y, '+b');hold on
% plot([ill_r(1) x], [ill_r(2) y], 'm-');hold on
% 
% axis([0 2 0 2])