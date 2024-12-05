function CI = Color_Index5( ill1, ill2, ill_r1, ill_r2 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define the line from ill1 to ill2 L1
m1 = (ill1(2) - ill2(2)) / (ill1(1) - ill2(1));
b1 = ill1(2) - m1*ill1(1);

%% Define the perpendicular line that passes through ill_r1
m2 = -1/m1;
b2 = ill_r1(2) - m2*ill_r1(1);

%% Find the intersection L1 & L2
x1 = (b2 -b1) / (m1-m2);
y1 = m1*x1 + b1;

%% Define the perpendicular line that passes through ill_r2
m2 = -1/m1;
b2 = ill_r2(2) - m2*ill_r2(1);

%% Find the intersection L1 & L2
x2 = (b2 -b1) / (m1-m2);
y2 = m1*x2 + b1;

%% Compute the CI metric
CI = norm( [x1,y1] - [x2,y2] ) / norm( ill1 - ill2 );

disp 'Negative condition disable'
% if x2 < x1
%     CI = 0;
% end

% figure;
% plot([ill1(1) ill2(1)], [ill1(2) ill2(2)], 'co-');hold on
% plot(ill_r(1), ill_r(2), '*k');hold on
% plot(x, y, '+b');hold on
% plot([ill_r(1) x], [ill_r(2) y], 'm-');hold on
% 
% axis([0 2 0 2])