%read_CYGNSS
clear; 

files = dir('/data/rlohman/HyperA/CYGNSS/level3/2018/*/*nc');
nd    =length(files);
dn0   =datenum(2018,01,01);
load('/data/pmb229/other/clearcuttingTStest/WorldHiVectors.mat'); 



for i=1:nd
    name=[files(i).folder '/' files(i).name];
    doy=name(end-5:end-3);
    data(i).dn=dn0+str2num(doy)-1;
    data(i).lon=ncread(name,'longitude');
    data(i).lat=ncread(name,'latitude');
    data(i).data=ncread(name,'SM_daily');
end




%     for i=200:220
%         pcolor(data(i).lon,data(i).lat,data(i).data); hold on; 
%         plot(lon, lat, 'k'); hold off; 
%         shading flat
%         axis([20 60 10 30])
%         title(datestr(data(i).dn))
%         Mov(i)=getframe;
%     end

%% time series
p = [55.105, 19.315]; 
%p = [54.259, 19.116]; 
va = []; 
idxa = [];
for i = 140:170
    ln = data(i).lon; 
    lt = data(i).lat; 
    da = data(i).data; 
    da = da(:); 
    [idx,val] = dsearchn([ln(:) lt(:)], p); 
    va = [va; da(idx)]; 
    idxa = [idxa; idx];
end
dts = 737061:737061+364;
dtsi = dts(140:170);

ix = ~isnan(va); 
va = va(ix); 
dtsi = dtsi(ix); 

figure; 
plot(dtsi, va, '.-'); 
datetick;
close; 

%%
cd /data/pmb229/other/SM/cygnss/
clearvars -except data
dd = datenum('2018/07/23')-737060;
lat2 = double(data(dd).lat);
lon2 = double(data(dd).lon);
d    = double(data(dd).data);

F = scatteredInterpolant(lon2(:), lat2(:), d(:), 'nearest');
[xv yv] = meshgrid(min(min(lon2)):0.2:max(max(lon2)), min(min(lat2)):0.2:max(max(lat2)));
vq = F(xv(:), yv(:));
vv = reshape(vq, size(xv,1), size(xv,2));


% close all; 
% figure; hold on;
% subplot(2,1,1); hold on; 
% plot(lon, lat, 'k');
% pcolor(xv, yv, vv); %shading flat;
% subplot(2,1,2); hold on; 
% plot(lon, lat, 'k');
% pcolor(lon2, lat2, d); %shading flat;
% colormap jet
% %axis([20 60 10 30]);
% 
% linkaxes; 
% 
% axis([50 60 15 20]);


% geotiff
d99 = vv; 
d99(isnan(d99)) = -99;
lonlim = [min(xv(:)) max(xv(:))]; 
latlim = [min(yv(:)) max(yv(:))]; 
s = size(d99);
R = georefcells(latlim, lonlim, s);
%geotiffwrite('cygnss_7-23-2018.tif', d99, R);
%system('mv /data/pmb229/other/SM/cygnss_* /data/pmb229/other/SM/cygnss/');

































