
im = double(imread(['D:\Matlab\CC2\rsrc\results_ColorConstancy2\' ...
    'screenshots\f_ILI1.png']))./255;
iwp = double(imread(['D:\Matlab\CC2\rsrc\results_ColorConstancy2\' ...
    'screenshots\f_wp_ILI1.png']))./255;

regoi_wp = reshape(iwp(266:373, 296:422, :), [], 3);

aux = 0:4;
regoi = [];
r0 = 351;
r1 = 377;
c0 = 114;
c1 = 145;
step = 124;

% regoi_wp = reshape(iwp(299:398, 312:417, :), [], 3);
% 
% aux = 0:4;
% regoi = [];
% r0 = 370;
% r1 = 392;
% c0 = 166;
% c1 = 196;
% step = 105;

for i=1:length(aux)
    regoi = [regoi;r0 r1 c0+aux(i)*step c1+aux(i)*step];
end

regoi_stat = [];
mask = zeros(size(im));
for i=1:length(aux)
    aux_im = ...
        reshape(im(regoi(i, 1):regoi(i, 2), ...
        regoi(i, 3):regoi(i, 4),:), [], 3);
    mask(regoi(i, 1):regoi(i, 2), ...
        regoi(i, 3):regoi(i, 4),:) = 1;
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
colormap default
for ii = 1:length(aux)
    for jj = 1:length(aux)
        text(ii, jj, num2str(round(deltaE_kulas(jj,ii),2)), ...
            'FontSize', 18);
    end
end
axis off