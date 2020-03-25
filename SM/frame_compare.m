% T28 info
info1 = importdata('/data/rlohman/Sentinel/Saudi/T28/geo_VV/rel_20170909_4r_4a.cor.geo.vrt'); 
fid   = fopen('/data/rlohman/Sentinel/Saudi/results_TS_VV_T28/20180524.time0'); 
fid2  = fopen('/data/rlohman/Sentinel/Saudi/results_TS_VV_T28/20180524.mag0'); 
fid3  = fopen('/data/rlohman/Sentinel/Saudi/results_TS_VV_T28/20180524.timeerr');

nx28  = 3042; 
ny28  = 4737; 
t28   = fread(fid, [nx28 ny28], 'real*4')'; fclose(fid);
m28   = fread(fid2, [nx28 ny28], 'real*4')'; fclose(fid2);
terr28   = fread(fid3, [nx28 ny28], 'real*4')'; fclose(fid3);
terr28(terr28>15) = NaN;
t28(isnan(terr28)) = NaN;
% directory
d     = dir('/data/rlohman/Sentinel/Saudi/T28/geo_VV/rel_*4r_4a.cor.geo');
d28   = {d.name};
% get dates
for i = 1:length(d28)
        n1    = cell2mat(d28(i));
        dn28(i) = datenum(n1(5:12), 'yyyymmdd');
end

% T130 info
info2 = importdata('/data/rlohman/Sentinel/Saudi/T130/geo_VV/rel_20170904_4r_4a.cor.geo.vrt'); 
fid   = fopen('/data/rlohman/Sentinel/Saudi/results_TS_VV_T130/20180524.time0'); 
fid2  = fopen('/data/rlohman/Sentinel/Saudi/results_TS_VV_T130/20180524.mag0');
fid3  = fopen('/data/rlohman/Sentinel/Saudi/results_TS_VV_T130/20180524.timeerr');
nx130 = 3159; 
ny130 = 3762; 
t130  = fread(fid, [nx130 ny130], 'real*4')'; fclose(fid);
m130  = fread(fid2, [nx130 ny130], 'real*4')'; fclose(fid2);
terr130  = fread(fid3, [nx130 ny130], 'real*4')'; fclose(fid3);
terr130(terr130>15) = NaN;
t130(isnan(terr130)) = NaN; 

% directory
d     = dir('/data/rlohman/Sentinel/Saudi/T130/geo_VV/rel_*4r_4a.cor.geo');
d130  = {d.name};
% get dates
for i = 1:length(d130)
        n1    = cell2mat(d130(i));
        dn130(i) = datenum(n1(5:12), 'yyyymmdd');
end

% Overlap
info3 = importdata('/data/rlohman/Sentinel/Saudi/T130_T28_resamp/rel_20170921.cor.geo.vrt');
fid   = fopen('/data/rlohman/Sentinel/Saudi/results_TS_T130_T28_VV/20180524.time0');
fid2  = fopen('/data/rlohman/Sentinel/Saudi/results_TS_T130_T28_VV/20180524.mag0');
fid3  = fopen('/data/rlohman/Sentinel/Saudi/results_TS_T130_T28_VV/20180524.timeerr');
nx    = 822;
ny    = 2538;
t0    = fread(fid, [nx ny], 'real*4')'; fclose(fid);
m0    = fread(fid2, [nx ny], 'real*4')'; fclose(fid2);
terr0 = fread(fid3, [nx ny], 'real*4')'; fclose(fid3);
% directory
d     = dir('/data/rlohman/Sentinel/Saudi/T130_T28_resamp/rel_*.cor.geo');
d     = {d.name};
% get dates
for i = 1:length(d)
        n1    = cell2mat(d(i));
        dn(i) = datenum(n1(5:12), 'yyyymmdd');
end

% storm date
sd28  = find(dn28 == datenum('20180531', 'yyyymmdd'));
sd130 = find(dn130 == datenum('20180526', 'yyyymmdd'));
sd    = find(dn == datenum('20180526', 'yyyymmdd'));


% import lat/long info
%T28
l3     = cell2mat(info1(3)); 
af     = strfind(l3, '>'); af2 = strfind(l3, '<'); cf = strfind(l3, ','); 
lonmin = str2double(l3(af(1)+1:cf(1)-1)); lond   = str2double(l3(cf(1)+1:cf(2)-1)); 
latmin = str2double(l3(cf(3)+1:cf(4)-1)); latd   = str2double(l3(cf(5)+1:af2(2)-1)); 

lonv28   = [lonmin:lond:lonmin+((nx28-1)*lond)]'; 
latv28   = [latmin:latd:latmin+((ny28-1)*latd)]'; 

%T130
l3     = cell2mat(info2(3)); 
af     = strfind(l3, '>'); af2 = strfind(l3, '<'); cf = strfind(l3, ',');
lonmin = str2double(l3(af(1)+1:cf(1)-1)); lond   = str2double(l3(cf(1)+1:cf(2)-1));
latmin = str2double(l3(cf(3)+1:cf(4)-1)); latd   = str2double(l3(cf(5)+1:af2(2)-1));
 
lonv130   = [lonmin:lond:lonmin+((nx130-1)*lond)]'; 
latv130   = [latmin:latd:latmin+((ny130-1)*latd)]'; 

%overlap
l3     = cell2mat(info3(3)); 
af     = strfind(l3, '>'); af2 = strfind(l3, '<'); cf = strfind(l3, ',');
lonmin = str2double(l3(af(1)+1:cf(1)-1)); lond   = str2double(l3(cf(1)+1:cf(2)-1));
latmin = str2double(l3(cf(3)+1:cf(4)-1)); latd   = str2double(l3(cf(5)+1:af2(2)-1));

lonv   = [lonmin:lond:lonmin+((nx-1)*lond)]';
latv   = [latmin:latd:latmin+((ny-1)*latd)]';

% point of interest
%xpt = 55.1208; 
%ypt = 19.347; 
xpt = 55.632; 
ypt = 18.192; 
xpt = 55.172; 
ypt = 19.512; 
x    = [1: 0.1:100]; 


[idx28x,b]  = dsearchn(lonv28, xpt);
[idx28y,b]  = dsearchn(latv28, ypt);
[idx130x,b] = dsearchn(lonv130, xpt);
[idx130y,b] = dsearchn(latv130, ypt);
[idxx,b]    = dsearchn(lonv, xpt);
[idxy,b]    = dsearchn(latv, ypt);

 % T28 
 % pull out rel cor values at pixels of interest
cd /data/rlohman/Sentinel/Saudi/T28/geo_VV/
for j = 1:length(xpt);
for i = 1:length(d28)
        ni = cell2mat(d28(i));
        fid = fopen(ni, 'r');
        fseek(fid, (nx28*(idx28y(j))+idx28x(j))*4, -1);
        p = fread(fid, 1, 'real*4');
        fclose(fid);
        v28(i,j) = p;
end
end
a28  = -log(diag(m28(idx28y+1, idx28x+1))); 
b28  = 1./diag(t28(idx28y+1, idx28x+1));
f28  = a28*exp(-b28*x);



% T130
% pull out rel cor values at pixels of interest
 cd /data/rlohman/Sentinel/Saudi/T130/geo_VV/
 for j = 1:length(xpt);
 for i = 1:length(d130)
         ni = cell2mat(d130(i));
         fid = fopen(ni, 'r');
         fseek(fid, (nx130*(idx130y(j))+idx130x(j))*4, -1);
         p = fread(fid, 1, 'real*4');
         fclose(fid);
         v130(i,j) = p;
 end
 end
a130  = -log(diag(m130(idx130y+1, idx130x+1))); 
b130  = 1./diag(t130(idx130y+1, idx130x+1));
f130  = a130*exp(-b130*x);


% Overlap
% pull out rel cor values at pixels of interest
cd /data/rlohman/Sentinel/Saudi/T130_T28_resamp/
 for j = 1:length(xpt);
 for i = 1:length(d)
         ni = cell2mat(d(i));
         fid = fopen(ni, 'r');
         fseek(fid, (nx*(idxy(j))+idxx(j))*4, -1);
         p = fread(fid, 1, 'real*4');
	 fclose(fid);
         v(i,j) = p;
 end
 end
a0  = -log(diag(m0(idxy+1, idxx+1)));
b0  = 1./diag(t0(idxy+1, idxx+1));
f0  = a0*exp(-b0*x);



%%

 cm = [0.8 0.5 0.1; 0.8 0.1 0.1; 0.2 0.2 0.8];
dstrm = datenum('20180524', 'yyyymmdd'); 

%close all; 
figure('Name', [num2str(xpt) ', ' num2str(ypt)]); hold on; box on; 
plot(dn28(sd28-3:sd28+10), v28(sd28-3:sd28+10), '.', 'markersize', 15, 'color', cm(1,:));
plot(dn130(sd130-3:sd130+10), v130(sd130-3:sd130+10), '.', 'markersize', 15, 'color', cm(2,:));
plot(dn(sd-5:sd+20), v(sd-5:sd+20), '.', 'markersize', 15, 'color', cm(3,:));

h1 = plot(x+dstrm-1, 1-f28, '-', 'linewidth', 2, 'color', cm(1,:));
h2 = plot(x+dstrm-1, 1-f130, '-', 'linewidth', 2, 'color', cm(2,:));
h3 = plot(x+dstrm-1, 1-f0, '-', 'linewidth', 2, 'color', cm(3,:));
 
legend([h1, h2, h3], 'T28', 'T130', 'Overlap', 'location', 'southeast');

%keyboard;
close; 



%
[latg28, long28]   = meshgrid(latv28, lonv28);
[latg130, long130] = meshgrid(latv130, lonv130); 
[latg, long]       = meshgrid(latv, lonv);

k28   = dsearchn([long28(:) latg28(:)], [lonmin latmin]); 
k130  = dsearchn([long130(:) latg130(:)], [lonmin latmin]); 
dd28  = [long28(k28) latg28(k28)]-[lonmin latmin]; % plot(long(:)+dd(1), latg(:)+dd(2), 'o'); 
dd130 = [long130(k130) latg130(k130)]-[lonmin latmin];

[i28,j28]   = ind2sub(size(long28), k28);
[i130,j130] = ind2sub(size(long130), k130);

pm1 = 1;
pm2 = -2; 
long28n  = long28(i28-pm1:i28+nx+pm2, j28-pm1:j28+ny+pm2); 
latg28n  = latg28(i28-pm1:i28+nx+pm2, j28-pm1:j28+ny+pm2);
long130n = long130(i130-pm1:i130+nx+pm2, j130-pm1:j130+ny+pm2);
latg130n = latg130(i130-pm1:i130+nx+pm2, j130-pm1:j130+ny+pm2);
%t28n     = t28(j28-1:j28+ny-2, i28-1:i28+nx-2); 
%t130n    = t130(j130-1:j130+ny-2, i130-1:i130+nx-2); 
t28n     = t28(j28-pm1:j28+ny+pm2, i28-pm1:i28+nx+pm2)';
t130n    = t130(j130-pm1:j130+ny+pm2, i130-pm1:i130+nx+pm2)';

v28      = griddata(long28n, latg28n, t28n, long, latg); 
v130     = griddata(long130n, latg130n, t130n, long, latg); 


%% map view figure comparison

terr0(terr0>10) = NaN; 

diff1 = t0-t28n';     
diff2 = t0-t130n'; 
diff3 = t28n'-t130n'; 
diff4 = t0-v28';      diff4(isnan(terr0)) = NaN; 
diff5 = t0-v130';     diff5(isnan(terr0)) = NaN; 
diff6 = v28'-v130';   diff6(isnan(terr0)) = NaN;

g = 0.8; 
clims = [-10 10]; 


%figure('units','normalized','outerposition',[0 .9 1 .5]); hold on; 
%subplot(1,3,1); hold on; 
%surf(lonv, latv, t0-t28n'); view(2); shading flat;
%caxis(clims); 
%colorbar
%title('Overlap-T28');
%set(gca,'Color',[1 1 1]*g)

%subplot(1,3,2); hold on; 
%surf(lonv, latv, t0-t130n'); view(2); shading flat;
%caxis(clims); 
%colorbar
%title('Overlap-T130'); 
%set(gca,'Color',[1 1 1]*g)

%subplot(1,3,3); hold on; 
%surf(lonv, latv, t28n'-t130n'); view(2); shading flat;
%caxis(clims); 
%title('T28-T130'); 
%colorbar; 
%set(gca,'Color',[1 1 1]*g)

%colormap bluewhitered
%dcmObj = datacursormode;
%set(dcmObj, 'UpdateFcn', @datacursorprecision);
%linkaxes; 

figure('units','normalized','outerposition',[0 .6 1 .5]); hold on;
subplot(1,3,1); hold on;
surf(lonv, latv, diff4); view(2); shading flat;
caxis(clims);
colorbar
title('Overlap-T28');
set(gca,'Color',[1 1 1]*g)

subplot(1,3,2); hold on; 
surf(lonv, latv, diff5);view(2); shading flat;
caxis(clims);
colorbar
title('Overlap-T130');
set(gca,'Color',[1 1 1]*g)

subplot(1,3,3); hold on; 
surf(lonv, latv, diff6);view(2); shading flat;
caxis(clims);
title('T28-T130');
colorbar;
set(gca,'Color',[1 1 1]*g)

colormap bluewhitered
dcmObj = datacursormode;
set(dcmObj, 'UpdateFcn', @datacursorprecision);
linkaxes; 
axis([55 55.7 17.7 19.8]);




















