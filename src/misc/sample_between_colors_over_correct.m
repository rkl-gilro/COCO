function [c_rgb, c_lab] = sample_between_colors_over_correct(c1, c2, n, n_o)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Computes competitors between c1 and c2
%% n is the number of samples between c1 and c2 (both included)
%% n_o is the number of overconstant values 
%% c1 is always the 'neutral' value
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norm(c1-c2)
deltaE00(c1', c2')

t = c2-c1;
step = 1/(n-1);

for i=1:n_o
    
    c_lab(i, :) = c1 - (i)*1.25*step*t;
    c_rgb(i, :) = lab2rgb(c_lab(i, :), 'WhitePoint', [1 1 1], 'ColorSpace', 'linear-rgb');
    
end

for i=1:n
    
    c_lab(n_o+i, :) = c1 + (i-1)*step*t;
    c_rgb(n_o+i, :) = lab2rgb(c_lab(n_o+i, :), 'WhitePoint', [1 1 1], 'ColorSpace', 'linear-rgb');
    disp(['Eab: ' , num2str(norm(c1-c_lab(n_o+i, :)))])
    disp(['E00: ' , num2str(deltaE00(c1', c_lab(n_o+i, :)'))])
end

