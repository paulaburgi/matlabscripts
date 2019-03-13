%% Create two axes
figure; hold on; 
ax1 = axes;
pcolor(ax1, mag(1e3:1.6e3, 1e3:1.6e3)); shading flat; 
caxis([0 1e7])
axis off
alpha 0.9
ax2 = axes;
pcolor(ax2, phs(1e3:1.6e3, 1e3:1.6e3)); shading flat; 
alpha 0.01

%% Link them together
linkaxes([ax1,ax2])
%% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%% Give each one its own colormap
colormap(ax1,'jet')
colormap(ax2,'gray')
%% Then add colorbars and get everything lined up
%set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
