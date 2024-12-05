
%% CC2 Illuminant experiment
im_cubes = double(imread('Data/LocalSurround/ScreenShot00012.png'))./255;
im_wp = double(imread('Data/LocalSurround/ScreenShot00009.png'))./255;

reg_wp = [153 265;608 759];

for i=1:3
    wp_value(i) = mean2...
    (im_wp(reg_wp(1,1):reg_wp(1,2), reg_wp(2,1):reg_wp(2,2), i));
end

reg_cubes = [308 342;260 306];
step = 208;
for j=1:5
    for i=1:3
        cube_value(j, i) = mean2...
            (im_cubes(reg_cubes(1,1):reg_cubes(1,2), ...
            reg_cubes(2,1):reg_cubes(2,2)+(j-1)*step, i));
    end
end