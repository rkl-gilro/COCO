function plot_densitymaps(xyz)

%% Compute density
% Put points into 3D bins; xyzBinNum is an nx3 matrix containing
% the bin ID for n values in xyz for the [x,y,z] axes.
nBins = 8;  % number of bins
sz = 100;
xbins = linspace(min(xyz(:,1)),max(xyz(:,1))*1,nBins+1);
ybins = linspace(min(xyz(:,2)),max(xyz(:,2))*1,nBins+1);
zbins = linspace(min(xyz(:,3)),max(xyz(:,3))*1,nBins+1);
xyzBinNum = [...
    discretize(xyz(:,1),xbins), ...
    discretize(xyz(:,2),ybins), ...
    discretize(xyz(:,3),zbins), ...
    ];
% bin3D is a mx3 matrix of m unique 3D bins that appear 
% in xyzBinNum, sorted.  binNum is a nx1 vector of bin
% numbers identifying the bin for each xyz point. For example,
% b=xyz(j,:) belongs to bins3D(b,:).
[bins3D, ~, binNum] = unique(xyzBinNum, 'rows');
% density is a mx1 vector of integers showing the number of 
% xyz points in each of the bins3D. To see the number of points
% in bins3D(k,:), density(k).  
density = histcounts(binNum,[1:size(bins3D,1),inf])'; 
% Compute bin centers
xbinCnt = xbins(2:end)-diff(xbins)/2;
ybinCnt = ybins(2:end)-diff(ybins)/2;
zbinCnt = zbins(2:end)-diff(zbins)/2;
%% Plot raw data
fig = figure();
tiledlayout('flow')
nexttile()
scatter3(xyz(:,1), ...
    xyz(:,2), ...
    xyz(:,3),sz)
title('Raw data')
grid on
box on
%% Plot bubblechart3
nexttile()
bubblechart3(...
    xbinCnt(bins3D(:,1)), ...
    ybinCnt(bins3D(:,2)), ...
    zbinCnt(bins3D(:,3)), ...
    density*50, ...
    density, ...
    'MarkerFaceAlpha', .3, ...
    'MarkerEdgeAlpha', .3)
title('bubblechart3')
%% Plot scatter3
nexttile()
scatter3(...
    xbinCnt(bins3D(:,1)), ...
    ybinCnt(bins3D(:,2)), ...
    zbinCnt(bins3D(:,3)), ...
    density*5+5, ...
    density, 'filled', ...
    'MarkerFaceAlpha',.6)
title('scatter3')
box on
%% Plot marginal histograms
% Histograms are normalized and are not related to 
% the axis ticks though their relative height within
% the axes shows the population density (ie, 30% of 
% the axis hight means that the bin contains 30% of
% the data.  
nexttile()
scatter3(...
    xyz(:,1), ...
    xyz(:,2), ...
    xyz(:,3), sz)
grid on
box on
hold on
xcount = histcounts(xyz(:,1),xbins)/size(xyz,1)*range(zlim())+min(zlim()); 
ycount = histcounts(xyz(:,2),ybins)/size(xyz,1)*range(zlim())+min(zlim()); 
% x density plotted on xz plane
patch(repelem(xbins(:),2,1), ...
    repmat(max(ybins),numel(xbins)*2,1), ...
    [min(zbins);repelem(xcount(:),2,1);min(zbins)], ...
    'r','FaceAlpha',.25)
% y density plotted on yz plane
patch(repmat(max(xbins),numel(xbins)*2,1), ...
    repelem(ybins(:),2,1), ...
    [min(zbins);repelem(ycount(:),2,1);min(zbins)], ...
    'r','FaceAlpha',.25)
title('Marginal histograms')
%% Plot marginal surface plots
% All 3 surface planes use the same color limit. 
nexttile()
scatter3(...
    xyz(:,1), ...
    xyz(:,2), ...
    xyz(:,3), ...
    sz, 'w', 'filled', ...
    'MarkerFaceAlpha', .5)
hold on
grid off
xlim([min(xbins), max(xbins)])
ylim([min(ybins), max(ybins)])
zlim([min(zbins), max(zbins)])
% xz plane density
[xm,zm] = ndgrid(xbins, zbins);
xzCount = histcounts2(xyz(:,1),xyz(:,3), xbins,zbins);
surf(xm,max(ylim())*ones(size(xm)),zm,xzCount,'FaceAlpha',.8)
% yz plane density
[ym,zm] = ndgrid(ybins, zbins);
yzCount = histcounts2(xyz(:,2),xyz(:,3), ybins,zbins);
surf(max(xlim())*ones(size(xm)),ym,zm,yzCount,'FaceAlpha',.8)
% xy plane density
[xm,ym] = ndgrid(xbins, ybins);
xyCount = histcounts2(xyz(:,1),xyz(:,2), xbins,ybins);
surf(xm,ym,min(zlim())*ones(size(xm)),xyCount,'FaceAlpha',.8)
set(gca,'LineWidth',3)
box on
maxCount = max([xyCount(:);xyCount(:);xzCount(:)]);
set(gca,'colormap',gray(256),'CLim',[0,maxCount])
cb = colorbar(); 
ylabel(cb,'Number of dots')
title('Density Planes')
%% Equate all axis limits
fig.UserData.hlink = linkprop(findobj(fig,'type','axes'),{'xlim','ylim','zlim'});