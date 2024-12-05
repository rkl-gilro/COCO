function XYZ = applyCalibration_rgb2xyz(file_calib, values)

%% Convert xyY to Unreal Standard
load(file_calib)                                                           % load calibration monXYZ & Gammas

x = (0:5:255)./255;
for ch = 1:3
    RGB_nl(:, ch) = interp1(x, radiometric(:, ch), values(:, ch));
end

XYZ = RGB_nl * monXYZ;


