function n_ill = name_illuminant( Illuminant )

experiments_ill = [0.2912 0.2705 0.2614; %neutral
    0.38 0.2441 0.1305;  %yellow
    0.233 0.3 0.2003;  %green
    0.2316 0.2833 0.4579;    %blue
    0.4152 0.2114 0.3961]; %red

diff = sum(experiments_ill - Illuminant, 2);

[~, ind] = min(abs(diff));

switch ind
    case 1
        n_ill = 'neutral';
        
    case 2
        n_ill = 'yellow';
        
    case 3
        
       n_ill = 'green';
       
    case 4
        n_ill = 'blue';
        
    case 5
        n_ill = 'red';
        
end