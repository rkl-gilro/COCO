% clc
% clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Main function for Matlab/Unreal connection measure primaries
%% D:\VR_Projects\CalibrationHMD unreal
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_filename = 'Calib_UnrealStandard521_Varjo_14_06_2024.mat';

N = 3;

pause(4);

% %% Create the connection
t = tcpclient('127.0.0.1', 8890);
fopen(t);


range = (0:5:255)./255;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Measure Red channel
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cont=1;
color=[];
i=1;
while i <= length(range)
    
    tic
    

    fwrite(t, "Value:" + range(i) + "," + 0 + "," + 0);
    a = fscanf(t, '%s\n');
    
    while ~strcmp(a, "SHOT")
        a = fscanf(t, '%s\n');
        fwrite(t, "Value:" + range(i) + "," + 0 + "," + 0);
    end
    
    pause(1);
    disp("Value:" + range(i) + "," + 0 + "," + 0)
    Red(i) = cs2000.measure;

    aux(cont) = Red(i);
    color = [color;Red(i).color.xyY'];
    
    cont = cont + 1;
    while(cont<N && i> 1)
        fwrite(t, "Value:" + range(i) + "," + 0 + "," + 0);
        a = fscanf(t, '%s\n');
        while ~strcmp(a, "SHOT")
            a = fscanf(t, '%s\n');
            fwrite(t, "Value:" + range(i) + "," + 0 + "," + 0);
        end
        aux(cont) = cs2000.measure;
        color = [color;aux(cont).color.xyY'];
        cont=cont+1;
    end
    if i>1
        [~, ind] = sort(sum(abs(color - median(color)),2));
        Red(i) = aux(ind(1));
    end
    disp(Red(i).color.xyY')
    
    
    t_time = toc;
    disp(['It took ', num2str(t_time), ' s']);
    disp '-------------------------------------------'
    i=i+1;
    cont = 1;
    color = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Measure Green channel
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cont=1;
color=[];
i=1;

while i <= length(range)
    
    tic
    
    fwrite(t, "Value:" + 0 + "," + range(i) +  "," + 0);
    a = fscanf(t, '%s\n');
    
    while ~strcmp(a, "SHOT")
        a = fscanf(t, '%s\n');
        fwrite(t, "Value:" + 0 + "," + range(i) +  "," + 0);
    end
    
    pause(1);
    disp("Value:" + 0 + "," + range(i) +  "," + 0)
    Green(i) = cs2000.measure;

    aux(cont) = Green(i);
    color = [color;Green(i).color.xyY'];
    
    cont = cont + 1;
    while(cont<N && i> 1)
        fwrite(t, "Value:" + 0 + "," + range(i) +  "," + 0);
        a = fscanf(t, '%s\n');
        while ~strcmp(a, "SHOT")
            a = fscanf(t, '%s\n');
            fwrite(t, "Value:" + 0 + "," + range(i) +  "," + 0);
        end
        aux(cont) = cs2000.measure;
        color = [color;aux(cont).color.xyY'];
        cont=cont+1;
    end
    if i>1
        [~, ind] = sort(sum(abs(color - median(color)),2));
        Green(i) = aux(ind(1));
    end
    disp(Green(i).color.xyY')

    
    t_time = toc;
    disp(['It took ', num2str(t_time), ' s']);
    disp '-------------------------------------------'
    i=i+1;
    cont = 1;
    color = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Measure Blue channel
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cont = 1;
color = [];
i = 1;
while i <= length(range)
    
    tic
    
    fwrite(t, "Value:" + 0 + "," + 0 + "," + range(i));
    a = fscanf(t, '%s\n');
    
    while ~strcmp(a, "SHOT")
        a = fscanf(t, '%s\n');
        fwrite(t, "Value:" + 0 + "," + 0 + "," + range(i));
    end
    
    pause(1);
    disp("Value:" + 0 + "," + 0 + "," + range(i))
    Blue(i) = cs2000.measure;

    aux(cont) = Blue(i);
    color = [color;Blue(i).color.xyY'];
    
    cont = cont + 1;
    while(cont<N && i> 1)
        disp("Value:" + 0 + "," + 0 + "," + range(i))
        a = fscanf(t, '%s\n');
        while ~strcmp(a, "SHOT")
            a = fscanf(t, '%s\n');
            disp("Value:" + 0 + "," + 0 + "," + range(i))
        end
        aux(cont) = cs2000.measure;
        color = [color;aux(cont).color.xyY'];
        cont=cont+1;
    end
    if i>1
        [~, ind] = sort(sum(abs(color - median(color)),2));
        Blue(i) = aux(ind(1));
    end
    disp(Blue(i).color.xyY')

    
    t_time = toc;
    disp(['It took ', num2str(t_time), ' s']);
    disp '-------------------------------------------'
    i=i+1;
    cont = 1;
    color = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Measure Gray channel
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cont = 1;
color = [];
i = 1;
while i <= length(range)
    
    tic
    
    fwrite(t, "Value:" + range(i) + "," + range(i) + "," + range(i));
    a = fscanf(t, '%s\n');
    
    while ~strcmp(a, "SHOT")
        a = fscanf(t, '%s\n');
        fwrite(t, "Value:" + range(i) + "," + range(i) + "," + range(i));
    end
    
    pause(1);
    disp("Value:" + range(i) + "," + range(i) + "," + range(i))
    Gray(i) = cs2000.measure;

    aux(cont) = Gray(i);
    color = [color;Gray(i).color.xyY'];
    
    cont = cont + 1;
    while(cont<N && i> 1)
        disp("Value:" + range(i) + "," + range(i) + "," + range(i))
        a = fscanf(t, '%s\n');
        while ~strcmp(a, "SHOT")
            a = fscanf(t, '%s\n');
            disp("Value:" + range(i) + "," + range(i) + "," + range(i))
        end
        aux(cont) = cs2000.measure;
        color = [color;aux(cont).color.xyY'];
        cont=cont+1;
    end
    if i>1
        [~, ind] = sort(sum(abs(color - median(color)),2));
        Gray(i) = aux(ind(1));
    end
    disp(Gray(i).color.xyY')

    
    t_time = toc;
    disp(['It took ', num2str(t_time), ' s']);
    disp '-------------------------------------------'
    i=i+1;
    cont = 1;
    color = [];
end
White = Gray(end);

values = range;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Validation predefined values
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load PredefinedRGB.mat
cont = 1;
color = [];
clear range

range = double(PredefinedRGB)./255;
for i = 1:size(PredefinedRGB, 1)
    
    tic
    
    fwrite(t, "Value:" + range(i, 1) + "," + range(i, 2) + "," ...
        + range(i, 3));
    a = fscanf(t, '%s\n');
    
    while ~strcmp(a, "SHOT")
        a = fscanf(t, '%s\n');
        fwrite(t, "Value:" + range(i, 1) + "," + range(i, 2) + ...
            "," + range(i, 3));
    end
    
    pause(1);
    disp("Value:" + range(i, 1) + "," + range(i, 2) + "," + range(i, 3))
    Validation_rand(i) = cs2000.measure;

    aux(cont) = Validation_rand(i);
    color = [color;Validation_rand(i).color.xyY'];
    
    cont = cont + 1;
    while(cont<N && i> 1)
        disp("Value:" + range(i, 1) + "," + range(i, 2) + ...
            "," + range(i, 3))
        a = fscanf(t, '%s\n');
        while ~strcmp(a, "SHOT")
            a = fscanf(t, '%s\n');
            disp("Value:" + range(i, 1) + "," + range(i, 2) + ...
                "," + range(i, 3))
        end
        aux(cont) = cs2000.measure;
        color = [color;aux(cont).color.xyY'];
        cont=cont+1;
    end
    if i>1
        [~, ind] = sort(sum(abs(color - median(color)),2));
        Validation_rand(i) = aux(ind(1));
    end
    disp(Validation_rand(i).color.xyY')
    cont=1;
    color = [];

    t_time = toc;
    disp(['It took ', num2str(t_time), ' s']);
    disp '-------------------------------------------'
    
end




save(save_filename, 'Red', 'Blue', 'Green', 'Gray', 'White', 'values',...
                    'Validation_rand', 'PredefinedRGB', '-v7.3');

                
fwrite(t,"DONE:0");