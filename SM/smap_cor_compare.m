%smap = geotiffread('/data/pmb229/other/SM/smapL2/times_2018_all2.tif'); 


% open timespan and error
%info = importdata('/data/rlohman/Sentinel/Saudi/T130/geo_VV/rel_20170904_4r_4a.cor.geo.vrt'); 
%fid  = fopen('/data/rlohman/Sentinel/Saudi/results_TS_VV_T130/20180524.time0'); 
%fid2 = fopen('/data/rlohman/Sentinel/Saudi/results_TS_VV_T130/20180524.timeerr');
%nx = 3159; 
%ny = 3762; 

%info = importdata('/data/rlohman/Sentinel/Saudi/T28/geo_VV/rel_20170909_4r_4a.cor.geo.vrt'); 
%fid  = fopen('/data/rlohman/Sentinel/Saudi/results_TS_VV/20180524.time0'); 
%fid2 = fopen('/data/rlohman/Sentinel/Saudi/results_TS_VV/20180524.timeerr'); 
%nx = 3042; 
%ny = 4737; 

info = importdata('/data/rlohman/Sentinel/Saudi/T130_T28_resamp/rel_20170921.cor.geo.vrt');
fid  = fopen('/data/rlohman/Sentinel/Saudi/results_TS_T130_T28_VV/20180524.time0');
fid2 = fopen('/data/rlohman/Sentinel/Saudi/results_TS_T130_T28_VV/20180524.timeerr');
nx = 822;
ny = 2538;

rmg  = fread(fid, [nx ny], 'real*4')'; 
rmg2 = fread(fid2, [nx ny], 'real*4')'; 
fclose(fid); 
fclose(fid2); 
idx = find(rmg2 > 50); 
rmg(idx) = NaN; 
rmg2(idx) = NaN; 
rmg(rmg == NaN) = -7; 

% import lat/long info
l3 = cell2mat(info(3)); 
af = strfind(l3, '>'); 
af2 = strfind(l3, '<'); 
cf = strfind(l3, ','); 
lonmin = str2double(l3(af(1)+1:cf(1)-1)); 
lond   = str2double(l3(cf(1)+1:cf(2)-1)); 
latmin = str2double(l3(cf(3)+1:cf(4)-1)); 
latd   = str2double(l3(cf(5)+1:af2(2)-1)); 

lonv   = [lonmin:lond:lonmin+((nx-1)*lond)]'; 
latv   = [latmin:latd:latmin+((ny-1)*latd)]'; 
[lonvg, latvg] = meshgrid(lonv, latv); 


% smap mags/times
load('mags_times2.mat');
long = s.long; 
latg = s.latg; 
ulat = s.ulat; % center of the cell. 
ulon = s.ulon; % center of the cell. 
times = s.times; 


%% histogram & map figure

% get all smap pixel boundaries
cboxall = []; 
for i = 2:length(ulat)-1
	for j = 2:length(ulon)-1
		clon = ulon(j);
		clat = ulat(i);
		% lat/lon of the limits of smap pixel
		lon1 = (clon-ulon(j-1))./2;
		lon2 = (clon-ulon(j+1))./2;
		lat1 = (clat-ulat(i+1))./2;
		lat2 = (clat-ulat(i-1))./2;
		cbox = [clat+lat2 clat+lat1 clon+lon2 clon+lon1];
		cboxall = [cboxall; cbox]; 
	end
end

% pixel 1 in smap to extract cor data
% ii = 14; 
% jj = 16; 
ii = 18; 
jj = 15; 
clon = ulon(jj); 
clat = ulat(ii); 
sm   = times(ii,jj);
% lat/lon of the limits of smap pixel
lon1 = (clon-ulon(jj-1))./2; 
lon2 = (clon-ulon(jj+1))./2; 
lat1 = (clat-ulat(ii+1))./2; 
lat2 = (clat-ulat(ii-1))./2; 
cbox = [clat+lat2 clat+lat1 clon+lon2 clon+lon1]; 

in = inpolygon(latvg, lonvg, cbox(1:2), cbox(3:4)); 
idx      = find(in == 1); 
[i1,j1]  = ind2sub(size(rmg), idx);
inrmg    = rmg(min(i1):max(i1), min(j1):max(j1)); 

load('/data/pmb229/other/clearcuttingTStest/WorldHiVectors.mat'); 
sarbox       = [54.744, 20.887; 57.076, 21.171; 57.859, 17.723; 55.493, 17.336; 54.743, 20.887];
sarbox2       = [52.6, 21.3; 55.01, 21.65; 55.945, 16.975; 53.62, 16.46; 52.6, 21.3];
bbox         = [14.4 27 49.7 60.9]; % latmin latmax, lonmin, lonmax


close all
figure('units', 'normalized', 'outerposition', [.1 .7 .3 .9]); hold on; 
subplot('position', [0.1 0.05 0.8 0.3]); hold on; box on;
histogram(inrmg(:)); 
yl = ylim; 
plot([sm sm], [0 yl(2)], 'k--'); 
box on; 
xlim([0 60]); 
ylabel('count')
xlabel('time scale (days)'); 
title('Histogram of cor-derived timescales in a single smap pixel');
legend('InSAR Coherence', ['SMAP (' num2str(round(sm)) ' day(s))'])
%close

% mask area not in overlap region
t = [sarbox(1,:); sarbox2(2,:); sarbox2(3,:); sarbox(4,:); sarbox(1,:)];
[xg, yg] = meshgrid(lonv, latv); 
test = inpolygon(xg(:), yg(:), t(:,1), t(:,2));
test2 = reshape(test, length(latv), length(lonv));
smm = sum(sum(isnan(rmg))); 
rmg(isnan(rmg)) = -8.5; 
test3  = rmg.*test2; 
test3(test3 == 0) = -10; 
rmg = test3; 

subplot('position', [0.1 0.43 0.8 0.52]); hold on; box on;
imagesc(lonv, latv, rmg); caxis([0 100]); set(gca, 'ydir', 'default'); hold on; 
for i = 1:length(cboxall)
         plot([cboxall(i,3) cboxall(i,4) cboxall(i,4) cboxall(i,3) cboxall(i,3)], [cboxall(i,1) cboxall(i,1) cboxall(i,2) cboxall(i,2) cboxall(i,1)], 'color', [1 1 1]*0.3);
end
plot3(lon, lat, ones(length(lat),1)*100, 'color', 0*[1 1 1]); 
plot([bbox(3) bbox(3) bbox(4) bbox(4) bbox(3)], [bbox(1) bbox(2) bbox(2) bbox(1) bbox(1)], 'k'); 
plot3(sarbox(:,1), sarbox(:,2), ones(length(sarbox),1)*100, 'k', 'linewidth', 2); 
plot3(sarbox2(:,1), sarbox2(:,2), ones(length(sarbox2),1)*100, 'k', 'linewidth', 2); 
plot([cbox(3) cbox(4) cbox(4) cbox(3) cbox(3)], [cbox(1) cbox(1) cbox(2) cbox(2) cbox(1)], 'k', 'linewidth', 2)
px1 = 263; py1 = 1128; 
px2 = 309; py2 = 1092; 
  cm = [0.8 0.5 0.1; 0.8 0.1 0.1]; 
plot(lonv(px1), latv(py1), '.', 'markersize', 10, 'color', cm(1,:)); 
plot(lonv(px2), latv(py2), '.', 'markersize', 10, 'color', cm(2,:)); 
% ylim([sarbox(3,2) sarbox(2,2)]);
% xlim([sarbox(1,1) sarbox(3,1)]); 
ylim([sarbox2(3,2) sarbox2(2,2)]);
xlim([sarbox2(1,1) sarbox2(3,1)]); 
% caxis([0 200]);
caxis([-10 50]); 
cmap = colormap(parula);
cmap = [[1 1 1]; 0.8*[1 1 1]; cmap]; 
%cmap = cmap+0.3;
%cmap = cmap./max(max(cmap));
colormap(cmap)
colorbar; 
axis([min(lonv)-1 max(lonv)+1 min(latv)-0.4 max(latv)+0.4]);
daspect([1 (111.32*cosd(23))/111.32 .1]) %1 degree long = 111.32*cos(lat) 
dcmObj = datacursormode;
set(dcmObj, 'UpdateFcn', @datacursorprecision);

%%

% close all;
%  figure; box on; 
% surf(ulon, ulat, times); caxis([0 25]); set(gca, 'ydir', 'default'); hold on; view(2); 
% % %surf(lonv, latv, rmg*100); shading flat; 
% % imagesc(ulon, ulat, times); caxis([0 5]); set(gca, 'ydir', 'default'); hold on; 
% % %imagesc(lonv, latv, rmg); caxis([0 100]); set(gca, 'ydir', 'default');
% plot3(lon, lat, ones(length(lat),1)*100, 'color', 0.5*[1 1 1]); 
% plot([bbox(3) bbox(3) bbox(4) bbox(4) bbox(3)], [bbox(1) bbox(2) bbox(2) bbox(1) bbox(1)], 'k'); 
% plot3(sarbox(:,1), sarbox(:,2), ones(length(sarbox),1)*100, 'k'); 
% plot3(sarbox2(:,1), sarbox2(:,2), ones(length(sarbox2),1)*100, 'k'); 
% plot([cbox(3) cbox(4) cbox(4) cbox(3) cbox(3)], [cbox(1) cbox(1) cbox(2) cbox(2) cbox(1)], 'k', 'linewidth', 2)
% ylim(bbox(1:2));
% xlim(bbox(3:4)); 
% colorbar;
% title('SMAP timescales')






















