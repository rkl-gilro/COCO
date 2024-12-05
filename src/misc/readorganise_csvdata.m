function [Chosen_lizard_RGB, Chosen_lizard_num, ...
          Shown_lizard_RGB, Shown_lizard_num] = readorganise_csvdata...
          ( A, colors_competitors, experiments_pos  )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define chosen lizard for each ill: 1) vector -RGB value
%%                                    2) integer -competitor number-
%% Do the same for all 5 lizards shown at each trial in S1,2,3,4,5
Chosen_lizard_RGB = nan(size(A, 1),3);
Chosen_lizard_num = nan(size(A, 1),1);

for i = 1:length(experiments_pos)
    
    pos = find(strcmp(experiments_pos{i},A.Experiment));

    Chosen_lizard_RGB(pos, :) = [A.Red(pos) A.Green(pos) A.Blue(pos)];
    
    Shown_lizard_RGB.S1(pos, :) = A.S1(pos);
    Shown_lizard_RGB.S2(pos, :) = A.S2(pos);
    Shown_lizard_RGB.S3(pos, :) = A.S3(pos);
    Shown_lizard_RGB.S4(pos, :) = A.S4(pos);
    Shown_lizard_RGB.S5(pos, :) = A.S5(pos);

    aux = colors_competitors{i};
    for j=1:length(pos)
        
        [~, ind] = min(sum( abs(Chosen_lizard_RGB(pos(j), :) - aux), 2));

        Chosen_lizard_num(pos(j)) = ind; 
        
        pos1 = str2num(cell2mat(Shown_lizard_RGB.S1(pos(j), :)));
        pos2 = str2num(cell2mat(Shown_lizard_RGB.S2(pos(j), :)));
        pos3 = str2num(cell2mat(Shown_lizard_RGB.S3(pos(j), :)));
        pos4 = str2num(cell2mat(Shown_lizard_RGB.S4(pos(j), :)));
        pos5 = str2num(cell2mat(Shown_lizard_RGB.S5(pos(j), :)));

        [~, ind1] = min(sum( abs(pos1 - aux), 2));
        [~, ind2] = min(sum( abs(pos2 - aux), 2));
        [~, ind3] = min(sum( abs(pos3 - aux), 2));
        [~, ind4] = min(sum( abs(pos4 - aux), 2));
        [~, ind5] = min(sum( abs(pos5 - aux), 2));

        Shown_lizard_num.S1(pos(j)) = ind1;
        Shown_lizard_num.S2(pos(j)) = ind2;
        Shown_lizard_num.S3(pos(j)) = ind3;
        Shown_lizard_num.S4(pos(j)) = ind4;
        Shown_lizard_num.S5(pos(j)) = ind5;

    end

end
