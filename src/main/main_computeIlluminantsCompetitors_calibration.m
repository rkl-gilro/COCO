
file_calib_unrealstandard = ...
    'Calib_UnrealStandard521_Varjo_13_06_2024_LUT_dE.mat';

Y = 40;                                                                    % value of Aston reference 50, in our experiments set to 30
samples = 4;                                                               % number of competitors -1

% neutral    = [0.3127, 0.3291, Y];
% blue       = [0.2614, 0.2702, Y];
% yellow     = [0.3823, 0.3837, Y];
achromatic = [0.3154, 0.3472, Y];

neutral = [0.31, 0.33, Y];
blue    = [0.25, 0.26, Y];
green   = [0.30,0.38, Y];
red     = [0.32,0.26, Y];
yellow  = [0.39, 0.39, Y];

illuminants = [neutral;blue;green;red;yellow;achromatic];

%% Load the white point of the calibration
xyz_wp = applyCalibration_rgb2xyz(file_calib_unrealstandard, [1 1 1])';

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

%% Define illuminants and gray in RGB engine
neutral_rgb = RGB_unreals(1, :);
blue_rgb = RGB_unreals(2, :);
yellow_rgb = RGB_unreals(3, :);

RefNeutral = RGB_unreals(4, :);

% ReflectedNoConstancy = RefNeutral .* IlluminantNeutral;

ReflectanceConstancy   = RefNeutral;
%% Define the Perfect constancy and 0 constancy in the neutral
xyz_aux = xyz_newvalidation(1, :) .* xyz_newvalidation(4, :);
xyz_perfectneutral = xyz_aux ./ xyz_newvalidation(1, 2);

xyz_reflectanceneutral = xyz_newvalidation(4, :);

% 0 constancy blue
xyz_aux = xyz_perfectneutral ./ xyz_newvalidation(2, :);
xyz_noconstancyblue = xyz_aux .* xyz_newvalidation(2, 2);

% 0 constancy yellow
xyz_aux = xyz_perfectneutral ./ xyz_newvalidation(3, :);
xyz_noconstancyyellow = xyz_aux .* xyz_newvalidation(3, 2);

lab_perfectneutral = xyz2lab(xyz_reflectanceneutral,...
            'WhitePoint', xyz_wp');
lab_zeroblue = xyz2lab(xyz_noconstancyblue,...
            'WhitePoint', xyz_wp');
lab_zeroyellow = xyz2lab(xyz_noconstancyyellow,...
            'WhitePoint', xyz_wp');


[~, aux ] = sample_between_colors_over_correct(lab_perfectneutral, ...
    lab_zeroblue, 3, 0);
[~, aux2] = sample_between_colors_over_correct(lab_perfectneutral, ...
    lab_zeroyellow, 3, 0);
all_colors_neutral_lab = [aux(end:-1:2, :); aux2];

aux_xyz = lab2xyz(all_colors_neutral_lab, 'WhitePoint',xyz_wp');
all_colors_neutral_rgb = applyCalibration_xyz2rgb...
    (file_calib_unrealstandard,aux_xyz);

