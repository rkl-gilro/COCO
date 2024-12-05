function lab_regoi = Labvalues_KulaILI(filename, wp_filename)

im = double(imread(filename))./255;
iwp = double(imread(wp_filename))./255;

regoi_wp = reshape(iwp(194:243, 349:432, :), [], 3);

aux = 0:4;
regoi = [];
r0 = 234;
r1 = 245;
c0 = 249;
c1 = 266;
step = 60;

for i=1:length(aux)
    regoi = [regoi;r0 r1 c0+aux(i)*step c1+aux(i)*step];
end

regoi_stat = [];
for i=1:length(aux)
    aux_im = ...
        reshape(im(regoi(i, 1):regoi(i, 2), ...
        regoi(i, 3):regoi(i, 4),:), [], 3);

    regoi_stat = ...
        [regoi_stat; ...
        mean(aux_im)];% ...
        % median(aux_im) ...
        % std(aux_im)];
end

lab_regoi = rgb2lab(regoi_stat, 'colorspace', 'linear-rgb',...
    'whitepoint', mean(regoi_wp));
