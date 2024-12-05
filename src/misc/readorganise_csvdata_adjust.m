function Chosen_lizard_RGB ...
     = readorganise_csvdata_adjust( A, experiments_pos )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define chosen lizard for each ill: 1) vector -RGB value
%%                                    2) integer -competitor number-
%% Do the same for all 5 lizards shown at each trial in S1,2,3,4,5
Chosen_lizard_RGB = nan(size(A, 1), 3);


for i = 1:length(experiments_pos)


    pos = find(strcmp(experiments_pos{i},A.Time));

    Chosen_lizard_RGB(pos, :) = table2array([A(pos, 7:9)]);


  
end


