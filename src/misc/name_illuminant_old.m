function n_ill = name_illuminant_old( Illuminant )

experiments_ill = [0.2467 0.2355 0.226; %neutral
    0.3248 0.2122 0.1117;  %yellow
    0.1891 0.2635 0.1698;  %green
    0.189 0.246 0.3998;    %blue
    0.3549 0.1819 0.3444]; %red

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