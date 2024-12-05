% Y = 30;                                                                    % 50-value of Aston reference, in our experiments set to ...
% neutral = [0.31, 0.33, Y];
% blue = [0.26, 0.26, Y];
% green = [0.30,0.38, Y]; 
% red = [0.33,0.26, Y]; 
% yellow = [0.40, 0.39, Y];
% 
% neutral_rgb = mtxyYtoRGB( neutral,  'PARSfull_newcalib2.mat')  
% blue_rgb = mtxyYtoRGB( blue,  'PARSfull_newcalib2.mat')  
% green_rgb = mtxyYtoRGB( green,  'PARSfull_newcalib2.mat')  
% red_rgb = mtxyYtoRGB( red,  'PARSfull_newcalib2.mat')  
% yellow_rgb = mtxyYtoRGB( yellow,  'PARSfull_newcalib2.mat')  

neutral_rgb = RGB_unreals(1, :);
blue_rgb = RGB_unreals(2, :);
green_rgb = RGB_unreals(3, :);
red_rgb = RGB_unreals(4, :);
yellow_rgb = RGB_unreals(5, :);
%% TODO compute possible candidates for local context from the diagonal
%% of the illuminantion axes 


neutral_lab = rgb2lab( neutral_rgb,'whitepoint', [1 1 1], 'colorspace', 'linear-rgb')  
blue_lab = rgb2lab( blue_rgb,'whitepoint', [1 1 1], 'colorspace', 'linear-rgb')  
green_lab = rgb2lab( green_rgb, 'whitepoint', [1 1 1], 'colorspace', 'linear-rgb' )  
red_lab = rgb2lab( red_rgb, 'whitepoint', [1 1 1], 'colorspace', 'linear-rgb')  
yellow_lab = rgb2lab( yellow_rgb, 'whitepoint', [1 1 1], 'colorspace', 'linear-rgb')  

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
% scatter3sph(X,Z,Y,'size', S, 'color', ([neutral_rgb;blue_rgb;green_rgb;red_rgb;yellow_rgb].*1.5).^.45, 'transp', 1)
light('Position',[1 1 1],'Style','infinit','Color',[1 1 1]);
lighting gouraud;
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
light('Position',[1 1 1],'Style','infinit','Color',[1 1 1]);
lighting gouraud;
set(gca,'FontSize',18,'LineWidth',2, 'FontName', 'LM Mono 12')

N = 4;

%% Convert illuminants to DKL
blue_candidate_dkl = mat_rgb2dkl * blue_candidates_rgb(N, :)';
green_candidate_dkl = mat_rgb2dkl * green_candidates_rgb(N, :)';
red_candidate_dkl = mat_rgb2dkl * red_candidates_rgb(N, :)';
yellow_candidate_dkl = mat_rgb2dkl * yellow_candidates_rgb(N, :)';

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
a_r(:, :, 1) = red_candidates_rgb(N, 1).*a_r(:,:, 1);
a_r(:, :, 2) = red_candidates_rgb(N, 2).*a_r(:, :,  2);
a_r(:, :, 3) = red_candidates_rgb(N, 3).*a_r(:, :,  3);

figure;imshow([a_g a_b a_y a_r])

