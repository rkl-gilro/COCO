function competitors = read_color_competitors( Illuminant )

%% The order of the illuminants Neutral, Blue, Yellow, Green and Red
   % [0.2467 0.2355 0.226;
   %  0.189 0.246 0.3998;
   %  0.3248 0.2122 0.1117;
   %  0.1891 0.2635 0.1698;
   %  0.3549 0.1819 0.3444];


switch Illuminant
    case 1
        competitors = ...
   [0.14633 0.211874 0.391033; ...
    0.151144 0.176076 0.250082; ...
    0.146317 0.145001 0.146794; ...
    0.220558 0.181979 0.141159; ...
    0.3095 0.224953 0.134474];
        
    case 3
        competitors = ...
    [0.14633 0.211874 0.391033; ...
    0.150717 0.187464 0.292596; ...
    0.1505 0.165213 0.211753; ...
    0.146317 0.145001 0.146794; ...
    0.14096 0.131106 0.10747];
        
    case 4
        
        competitors = ...
    [0.203881 0.138419 0.208684; ...
    0.18384 0.140917 0.186574; ...
    0.164656 0.143109 0.165961; ...
    0.146317 0.145001 0.146794; ...
    0.133108 0.146229 0.133336];
        
    case 2
        competitors = ...
    [0.3095 0.224953 0.134474; ...
    0.248507 0.195615 0.139052; ...
    0.194243 0.169009 0.14315; ...
    0.146317 0.145001 0.146794; ...
    0.114297 0.128621 0.149243];
        
    case 5
        competitors = ...
    [0.146323 0.27006 0.138576;  ...
    0.148351 0.223304 0.141886; ...
    0.148293 0.181712 0.144614; ...
    0.146317 0.145001 0.146794; ...
    0.143678 0.120498 0.148087];
        
end