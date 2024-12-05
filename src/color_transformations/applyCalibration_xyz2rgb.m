function RGB = applyCalibration_xyz2rgb(file_calib, values)

%% Convert xyY to Unreal Standard
load(file_calib)                                                           % load calibration monXYZ & Gammas

xyz_newvalidation = values;
invMatrix = inv(monXYZ);

RGB_nl = xyz_newvalidation * invMatrix;
x = (0:5:255)./255;
for ch = 1:3

    RGB(:, ch) = interp1(radiometric(:, ch), x, RGB_nl(:, ch));
    
end
