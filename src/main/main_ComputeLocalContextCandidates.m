normalize = 0;
file_calib_unrealstandard = ...
    'Calibration_UnrealStandard_HTCVive_14_03_2023_LUT_dE.mat';
Y = 30;                                                                    % value of Aston reference 50, in our experiments set to 30

neutral = [0.31, 0.33, Y];
blue    = [0.25, 0.26, Y];
green   = [0.30,0.38, Y];
red     = [0.32,0.26, Y];
yellow  = [0.39, 0.39, Y];
achromatic = [0.2991, 0.3243, Y];

illuminants = [neutral;blue;green;red;yellow;achromatic];

%% Convert xyY to Unreal Standard
load(file_calib_unrealstandard)                                             % load calibration monXYZ & Gammas

xyz_newvalidation = xyYToXYZ(illuminants')';
invMatrix = inv(monXYZ);

RGB_nl = xyz_newvalidation * invMatrix;
x = (0:5:255)./255;
for ch = 1:3

    RGB(:, ch) = interp1(radiometric(:, ch), x, RGB_nl(:, ch));
    
end
RGB_unreals = RGB;

neutral_rgb = RGB_unreals(1, :);
blue_rgb = RGB_unreals(2, :);
green_rgb = RGB_unreals(3, :);
red_rgb = RGB_unreals(4, :);
yellow_rgb = RGB_unreals(5, :);

RefNeutral = RGB_unreals(6, :);

neutral_lab = rgb2lab( neutral_rgb,'whitepoint', [1 1 1], ...
    'colorspace', 'linear-rgb')  
blue_lab = rgb2lab( blue_rgb,'whitepoint', [1 1 1], ...
    'colorspace', 'linear-rgb')  
green_lab = rgb2lab( green_rgb, 'whitepoint', [1 1 1], ...
    'colorspace', 'linear-rgb' )  
red_lab = rgb2lab( red_rgb, 'whitepoint', [1 1 1], ...
    'colorspace', 'linear-rgb')  
yellow_lab = rgb2lab( yellow_rgb, 'whitepoint', [1 1 1], ...
    'colorspace', 'linear-rgb')  

%% Compute diagonals
blue2neutral = neutral_lab - blue_lab;
green2neutral = neutral_lab - green_lab;
red2neutral = neutral_lab - red_lab;
yellow2neutral = neutral_lab - yellow_lab;

% compute angles
% green-blue
a_greenblue = acos(sum(blue2neutral(2:3).*green2neutral(2:3))./(norm(blue2neutral(2:3))*norm(green2neutral(2:3))));
rot_greenneutral =  [1 0 0;0 cos(a_greenblue/2) -sin(a_greenblue/2);0 sin(a_greenblue/2) cos(a_greenblue/2)]*green2neutral';

% green-yellow
a_greenyellow = acos(sum(yellow2neutral(2:3).*green2neutral(2:3))./(norm(yellow2neutral(2:3))*norm(green2neutral(2:3))));
rot_yellowneutral =  [1 0 0;0 cos(a_greenyellow/2) -sin(a_greenyellow/2);0 sin(a_greenyellow/2) cos(a_greenyellow/2)]*yellow2neutral';

% yellow-red
a_yellowred = acos(sum(yellow2neutral(2:3).*red2neutral(2:3))./(norm(yellow2neutral(2:3))*norm(red2neutral(2:3))));
rot_redneutral =  [1 0 0;0 cos(a_yellowred/2) -sin(a_yellowred/2);0 sin(a_yellowred/2) cos(a_yellowred/2)]*red2neutral';

% blue-red
a_bluered = acos(sum(blue2neutral(2:3).*red2neutral(2:3))./(norm(blue2neutral(2:3))*norm(red2neutral(2:3))));
rot_blueneutral =  [1 0 0;0 cos(a_bluered/2) -sin(a_bluered/2);0 sin(a_bluered/2) cos(a_bluered/2)]*blue2neutral';

%% Plot illuminants
S = 3;
figure;
X = [neutral_lab(1) blue_lab(1) green_lab(1) red_lab(1) yellow_lab(1)];
Y = [neutral_lab(2) blue_lab(2) green_lab(2) red_lab(2) yellow_lab(2)];
Z = [neutral_lab(3) blue_lab(3) green_lab(3) red_lab(3) yellow_lab(3)];
% scatter3sph(X,Z,Y,'size', S, 'color', ([neutral_rgb';blue_rgb';green_rgb';red_rgb';yellow_rgb'].*1.5).^.45, 'transp', 1)
% light('Position',[1 1 1],'Style','infinit','Color',[1 1 1]);
% lighting gouraud;
axis equal
hold on
plot3([neutral_lab(1) neutral_lab(1)], [neutral_lab(3) blue_lab(3)], [neutral_lab(2) blue_lab(2)], 'k-')
plot3([neutral_lab(1) neutral_lab(1)], [neutral_lab(3) green_lab(3)], [neutral_lab(2) green_lab(2)], 'k-')
plot3([neutral_lab(1) neutral_lab(1)], [neutral_lab(3) red_lab(3)], [neutral_lab(2) red_lab(2)], 'k-')
plot3([neutral_lab(1) neutral_lab(1)], [neutral_lab(3) yellow_lab(3)], [neutral_lab(2) yellow_lab(2)], 'k-')

plot3([neutral_lab(1) neutral_lab(1)], [neutral_lab(3) neutral_lab(3)-rot_greenneutral(3)], [neutral_lab(2) neutral_lab(2)-rot_greenneutral(2)], 'k--')
plot3([neutral_lab(1) neutral_lab(1)], [neutral_lab(3) neutral_lab(3)-rot_blueneutral(3)], [neutral_lab(2) neutral_lab(2)-rot_blueneutral(2)], 'k--')
plot3([neutral_lab(1) neutral_lab(1)], [neutral_lab(3) neutral_lab(3)-rot_yellowneutral(3)], [neutral_lab(2) neutral_lab(2)-rot_yellowneutral(2)], 'k--')
plot3([neutral_lab(1) neutral_lab(1)], [neutral_lab(3) neutral_lab(3)-rot_redneutral(3)], [neutral_lab(2) neutral_lab(2)-rot_redneutral(2)], 'k--')


% doublecheck!!
% a = [rot_yellowneutral(2);
% rot_yellowneutral(3)]
% 
% b = [neutral_lab(2)-yellow_lab(2);
%  neutral_lab(3)-yellow_lab(3)]
% 
% acos(sum(a.*b)/(norm(a)*norm(b)))


%% Compute candidates
t= 0:.25:1;

green_candidates0 = [neutral_lab(1) neutral_lab(1);
 neutral_lab(2) neutral_lab(2)-rot_greenneutral(2);
 neutral_lab(3) neutral_lab(3)-rot_greenneutral(3)];

green_candidates = green_candidates0(:, 1) + t.*(green_candidates0(:, 2) - green_candidates0(:, 1));
green_candidates_rgb = lab2rgb(green_candidates', 'whitepoint', [1 1 1], 'colorspace', 'linear-rgb');

blue_candidates0 = [neutral_lab(1) neutral_lab(1);
    neutral_lab(2) neutral_lab(2)-rot_blueneutral(2);
    neutral_lab(3) neutral_lab(3)-rot_blueneutral(3)];

blue_candidates = blue_candidates0(:, 1) + t.*(blue_candidates0(:, 2) - blue_candidates0(:, 1));
blue_candidates_rgb = lab2rgb(blue_candidates', 'whitepoint', [1 1 1], 'colorspace', 'linear-rgb');

yellow_candidates0 = [neutral_lab(1) neutral_lab(1);
    neutral_lab(2) neutral_lab(2)-rot_yellowneutral(2);
    neutral_lab(3) neutral_lab(3)-rot_yellowneutral(3)];

yellow_candidates = yellow_candidates0(:, 1) + t.*(yellow_candidates0(:, 2) - yellow_candidates0(:, 1));
yellow_candidates_rgb = lab2rgb(yellow_candidates', 'whitepoint', [1 1 1], 'colorspace', 'linear-rgb');

red_candidates0 = [neutral_lab(1) neutral_lab(1);
 neutral_lab(2) neutral_lab(2)-rot_redneutral(2);
 neutral_lab(3) neutral_lab(3)-rot_redneutral(3)];

red_candidates = red_candidates0(:, 1) + t.*(red_candidates0(:, 2) - red_candidates0(:, 1))
red_candidates_rgb = lab2rgb(red_candidates', 'whitepoint', [1 1 1], 'colorspace', 'linear-rgb');

X_c = [green_candidates(1, :) blue_candidates(1, :) yellow_candidates(1, :) red_candidates(1, :)];
Y_c = [green_candidates(2, :) blue_candidates(2, :) yellow_candidates(2, :) red_candidates(2, :)];
Z_c = [green_candidates(3, :) blue_candidates(3, :) yellow_candidates(3, :) red_candidates(3, :)];
% scatter3sph(X_c,Z_c,Y_c,'size', S, 'color', ([green_candidates_rgb;blue_candidates_rgb;yellow_candidates_rgb;red_candidates_rgb].^1.5).^.45, 'transp', 1)
% light('Position',[1 1 1],'Style','infinit','Color',[1 1 1]);
% lighting gouraud;
% set(gca,'FontSize',18,'LineWidth',2, 'FontName', 'LM Mono 12')


N = 4;
a_g = ones(200,200,3);
% a(1:100, :, 1) = green_candidates_rgb(2, 1).*a(1:100, :, 1);
% a(1:100, :, 2) = green_candidates_rgb(2, 2).*a(1:100, :, 2);
% a(1:100, :, 3) = green_candidates_rgb(2, 3).*a(1:100, :, 3);
a_g(:, :, 1) = green_candidates_rgb(N, 1).*a_g(:,:, 1);
a_g(:, :, 2) = green_candidates_rgb(N, 2).*a_g(:, :,  2);
a_g(:, :, 3) = green_candidates_rgb(N, 3).*a_g(:, :,  3);



a_b = ones(200,200,3);
% a(1:100, :, 1) = blue_candidates_rgb(2, 1).*a(1:100, :, 1);
% a(1:100, :, 2) = blue_candidates_rgb(2, 2).*a(1:100, :, 2);
% a(1:100, :, 3) = blue_candidates_rgb(2, 3).*a(1:100, :, 3);
a_b(:, :, 1) = blue_candidates_rgb(N, 1).*a_b(:,:, 1);
a_b(:, :, 2) = blue_candidates_rgb(N, 2).*a_b(:, :,  2);
a_b(:, :, 3) = blue_candidates_rgb(N, 3).*a_b(:, :,  3);


a_y = ones(200,200,3);
% a(1:100, :, 1) = yellow_candidates_rgb(2, 1).*a(1:100, :, 1);
% a(1:100, :, 2) = yellow_candidates_rgb(2, 2).*a(1:100, :, 2);
% a(1:100, :, 3) = yellow_candidates_rgb(2, 3).*a(1:100, :, 3);
a_y(:, :, 1) = yellow_candidates_rgb(N, 1).*a_y(:,:, 1);
a_y(:, :, 2) = yellow_candidates_rgb(N, 2).*a_y(:, :,  2);
a_y(:, :, 3) = yellow_candidates_rgb(N, 3).*a_y(:, :,  3);


a_r = ones(200,200,3);
% a(1:100, :, 1) = red_candidates_rgb(2, 1).*a(1:100, :, 1);
% a(1:100, :, 2) = red_candidates_rgb(2, 2).*a(1:100, :, 2);
% a(1:100, :, 3) = red_candidates_rgb(2, 3).*a(1:100, :, 3);
a_r(:, :, 1) = red_candidates_rgb(N, 1).*a_r(:, :, 1);
a_r(:, :, 2) = red_candidates_rgb(N, 2).*a_r(:, :, 2);
a_r(:, :, 3) = red_candidates_rgb(N, 3).*a_r(:, :, 3);

figure;imshow([a_g a_b a_y a_r])

%% Lab candiodates
lab_candidates = [green_candidates(:, N)';...
                  blue_candidates(:, N)';...
                  yellow_candidates(:, N)';...
                  red_candidates(:, N)'];
lab_illuminants = [neutral_lab;...
                  green_lab;...
                  blue_lab;...
                  yellow_lab;...
                  red_lab];

%% Convert to DKL
monxyY = XYZToxyY(monXYZ');
dkl_candidates = RGB2DKL([green_candidates_rgb(N, :);...
                  blue_candidates_rgb(N, :);...
                  yellow_candidates_rgb(N, :);...
                  red_candidates_rgb(N, :)], monxyY');
dkl_illuminants = RGB2DKL([neutral_rgb;...
                  green_rgb;...
                  blue_rgb;...
                  yellow_rgb;...
                  red_rgb], monxyY');
candidates_colors = [green_candidates_rgb(N, :);...
                  blue_candidates_rgb(N, :);...
                  yellow_candidates_rgb(N, :);...
                  red_candidates_rgb(N, :)];
illuminants_colors = [neutral_rgb;...
                  green_rgb;...
                  blue_rgb;...
                  yellow_rgb;...
                  red_rgb];
figure;scatter3(dkl_candidates(:, 1), dkl_candidates(:, 2), ...
    dkl_candidates(:, 3), 340, candidates_colors.^.45,'filled');hold on

plot3([dkl_illuminants(1, 1) dkl_illuminants(2, 1)], ...
    [dkl_illuminants(1, 2) dkl_illuminants(2, 2)], ...
    [dkl_illuminants(1, 3) dkl_illuminants(2, 3)],'k--');
plot3([dkl_illuminants(1, 1) dkl_illuminants(3, 1)], ...
    [dkl_illuminants(1, 2) dkl_illuminants(3, 2)], ...
    [dkl_illuminants(1, 3) dkl_illuminants(3, 3)],'k--');
plot3([dkl_illuminants(1, 1) dkl_illuminants(4, 1)], ...
    [dkl_illuminants(1, 2) dkl_illuminants(4, 2)], ...
    [dkl_illuminants(1, 3) dkl_illuminants(4, 3)],'k--');
plot3([dkl_illuminants(1, 1) dkl_illuminants(5, 1)], ...
    [dkl_illuminants(1, 2) dkl_illuminants(5, 2)], ...
    [dkl_illuminants(1, 3) dkl_illuminants(5, 3)],'k--');

scatter3(dkl_illuminants(:, 1), dkl_illuminants(:, 2), ...
    dkl_illuminants(:, 3), 340, illuminants_colors.^.45,'filled');hold on

figure;
[theta,rho,z] = cart2pol(dkl_candidates(:, 2), dkl_candidates(:, 3), ...
dkl_candidates(:, 1));
polarscatter(theta, rho, 540,candidates_colors.^.45,'filled');hold on
[theta,rho,z] = cart2pol(dkl_illuminants(:, 2), dkl_illuminants(:, 3), ...
dkl_illuminants(:, 1));
polarscatter(theta, rho, 1040,illuminants_colors.^.45,'filled')
rticklabels({})
thetaticklabels({});
text(pi*.5, .3, 'S-(L+M)', 'horiz', 'center', 'vert', 'top', ...
    'rotation', 0,'FontSize', 24, 'fontname','Times New Roman');
text(0, .3, 'L-M', 'horiz', 'center', 'rotation', 0,'FontSize', 24,...
    'fontname','Times New Roman');
set(gca,  'FontSize', 24, 'fontname','Times New Roman'); 
grid on
set(gcf,'renderer','Painters');


%% Lab plot
figure;scatter(lab_candidates(:, 2), ...
    lab_candidates(:, 3), 540, candidates_colors.^.45,'filled');hold on

scatter(lab_illuminants(:, 2), ...
    lab_illuminants(:, 3), 2040, illuminants_colors.^.45,'filled');hold on
grid off
axis off