%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Stablish spectroradiometer CS2000A connection and set it up
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cs2000 = CS2000;
% Synchronization
sync = CS2000_Sync('Internal', 90);
cs2000.setSync(sync);