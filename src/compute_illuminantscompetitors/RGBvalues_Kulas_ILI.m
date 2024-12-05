
im = double(imread('D:\Matlab\ColorCharacterization-main\rsrc\ILI3.png'))./255;
iwp = double(imread('D:\Matlab\ColorCharacterization-main\rsrc\wp_ILI3.png'))./255;

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

for i=1:length(aux)
    for j=i:length(aux)
        deltaE_kulas(i,j) = deltaE00(lab_regoi(i, :)', ...
                                     lab_regoi(j, :)');
        deltaE_kulas(j, i) = deltaE_kulas(i,j);
    end
end



figure;
imagesc(deltaE_kulas);
for ii = 1:length(aux)
    for jj = 1:length(aux)
        text(ii, jj, num2str(round(deltaE_kulas(jj,ii),2)), ...
            'FontSize', 18);
    end
end
axis off