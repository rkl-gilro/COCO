%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Compare same kulas through different ILI values
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = [1 2 6];

for i=1:length(k)
    filename = strcat(...
        'D:\Matlab\CC2\rsrc\cc2_screenshots\f_ILI',...
        num2str(k(i)),'.png');
    wp_filename = strcat(...
        'D:\Matlab\CC2\rsrc\cc2_screenshots\f_wp_ILI',...
        num2str(k(i)),'.png');
    lab_regoi(:, :, i) = Labvalues_KulaILI(filename, wp_filename);

end

for i=1:length(k)
    for j=i:length(k)
        deltaE_kulas(:,i,j) = deltaE00(lab_regoi(:, :, i)', ...
                                     lab_regoi(:, :, j)')';
        deltaE_kulas(:,j, i) = deltaE_kulas(:,i,j);
    end
end

figure;
color = colorcube;
marker = {'o', '*','+','<','diamond'};
for i=1:length(k)-1
plot(1:length(deltaE_kulas(:, 2, 1)), ...
    deltaE_kulas(:, i+1, 1), "Color", color(i*2, :), ...
    'Marker',marker{i},'LineWidth',2);
    hold on
text(length(deltaE_kulas(:, 2, 1))+.05,...
    deltaE_kulas(end, i+1, 1), strcat('ILI1--ILI', ...
    num2str(k(i+1))),...
    'Color',color(i*2, :),'FontSize',15);
end

axis([0 6 0 11])
xticks([1 2 3 4 5])
xticklabels({'K1','K2','K3','K4','K5'})

ylabel('DeltaE00');
set(gca,  'FontSize', 20, 'fontname','Times New Roman'); 
grid on
set(gcf,'renderer','Painters');

