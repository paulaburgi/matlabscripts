
%clear
close all

% Define directory and data location
%sumatra
pf_fol  = '/data/pmb229/isce/p446f7190_sumatra/'; 
cintfol = 'ints_SRTM/';  pol = 'HH'; 
datafol = [pf_fol 'data/']; 
intfol  = [pf_fol cintfol];
LSfol   = [pf_fol 'Landsat/']; 
im = sqrt(-1); 
    
intdirs = dir(intfol); 
intdirs = {intdirs.name};

% load int size, lat lon coords. 
    x  = importdata([intfol cell2mat(intdirs(3)) '/filt_topophase.flat.geo.vrt']);
    l1 = x{1}; 
    qf = strfind(l1, '"'); 
nx = str2double(l1(qf(1)+1:qf(2)-1)); 
ny = str2double(l1(qf(3)+1:qf(4)-1)); 
    l3 = x{3}; 
    aidx = strfind(l3, '>'); 
    aidx2 = strfind(l3, '<'); 
    cidx = strfind(l3, ','); 
ul_lon = str2double(l3(aidx(1)+1:cidx(1)-1)); 
dx_lon = str2double(l3(cidx(1)+1:cidx(2)-1));
ul_lat = str2double(l3(cidx(3)+1:cidx(4)-1));
dx_lat = str2double(l3(cidx(5)+1:aidx2(2)-1));
lr_lon = ul_lon+(dx_lon*nx); 
lr_lat = ul_lat+(dx_lat*ny);
ilon   = ul_lon:dx_lon:lr_lon-dx_lon; ilon = ilon'; 
ilat   = ul_lat:dx_lat:lr_lat-dx_lat; ilat = ilat'; 

% load Hansen Landsat tree loss (0-16, years after 2000)
ly       = [LSfol 'Hansen_GFC-2016-v1.4_lossyear_10N_100E.tif']; 
ly_saved = [ly(1:end-4) '_trimmed.mat']; 
if exist(ly_saved, 'file');
    load(ly_saved); 
else
    [h_all, lon, lat] = loadHansen(ly); 
    lonidx = find(lon > ul_lon & lon < lr_lon);   hlon = lon(lonidx)'; 
    latidx = find(lat < ul_lat & lat > lr_lat);   hlat = lat(latidx)'; 
    h      = h_all(latidx, lonidx); 
        save(ly_saved, 'h', 'hlon', 'hlat');
end

% find areas cleared in 2008 and 2009
idx08 = find(h == 8); 
idx09 = find(h == 9); 
idx0809 = [idx08; idx09]; 
idxn0809 = find(h ~= 8 & h ~= 9); 

h0809 = h; 
h0809(idxn0809) = NaN; 
h0809(idx08) = 1; 
h0809(idx09) = 2; 

hlonm   = repmat(hlon', length(hlat), 1); 
hlatm   = repmat(hlat, 1, length(hlon)); 
hlon0809           = hlonm; 
hlon0809(idxn0809) = NaN; 
hlon0809(idx08)    = 1; 
hlon0809(idx09)    = 2; 
hlat0809           = hlatm; 
hlat0809(idxn0809) = NaN; 
hlat0809(idx08)    = 1; 
hlat0809(idx09)    = 2; 

lon08 = hlonm(idx08); lon08 = lon08(:); 
lat08 = hlatm(idx08); lat08 = lat08(:); 
lon09 = hlonm(idx09); lon09 = lon09(:); 
lat09 = hlatm(idx09); lat09 = lat09(:); 

ilonm = repmat(ilon', length(ilat), 1); 
ilatm = repmat(ilat, 1, length(ilon)); 
ilonv = ilonm(:); 
ilatv = ilatm(:); 

% find the pixels of interferograms that correspond to the deforested
% pixels from Hansen data, TAKES FOREVER
% sidx08 = dsearchn([ilonv ilatv], [lon08 lat08]); % index of int file, with points
% sidx09 = dsearchn([ilonv ilatv], [lon09 lat09]); % index of int file, with points
%sidx_all = dsearchn([ilonv ilatv], [hlonm(:) hlatm(:)]); % index of int file, with points
%sidx_all_2 = dsearchn([hlonm(:) hlatm(:)], [ilonv ilatv]); % index of int file, with points
%save('dsearchn_results_sidx_all_mar05-2018', 'sidx_all', 'sidx_all_2'); 


% save('dsearchn_results_mar02-2018', 'sidx08', 'sidx09'); 

load([LSfol 'dsearchn_results_mar02-2018.mat']); 
load([LSfol 'dsearchn_results_sidx_all_mar05-2018']); 

i0809         = zeros(ny, nx);
i0809(sidx08) = 1; 
i0809(sidx09) = 2; 
i0809         = flipud(i0809); 




cmap = [0 0 1; 1 0 0; 1 1 1]; 
dx   = 1; 
figure; hold on; 
subplot(1,2,1); 
    imagesc(hlon(1:dx:end), hlat(1:dx:end), h0809(1:dx:end, 1:dx:end)); shading faceted; colormap(flipud(cmap)); 
subplot(1,2,2); 
    imagesc(ilon(1:dx:end), flipud(ilat(1:dx:end)), (i0809(1:dx:end, 1:dx:end))); shading faceted; colormap(flipud(cmap)); 
    linkaxes
    
filename   = ['/data/pmb229/isce/p446f7190_sumatra/ints_SRTM/' ...
              'int_090105_090220/filt_topophase.flat.geo']; 
fid         = fopen(filename, 'r','native');
[rmg,count] = fread(fid,[nx*2,ny],'real*4');
status      = fclose(fid);
real        = flipud((rmg(1:2:nx*2,1:ny))');
imag        = flipud((rmg(2:2:nx*2,1:ny))');
mag         = abs(real+im*imag);
phs         = angle(real+im*imag);

% figure; hold on; 
% ax1 = subplot(1,2,1); 
%     pcolor(ilon, ilat, i0809(1:1:end,1:1:end)); shading flat; colormap(ax1, flipud(cmap)); 
%     h=colorbar; 
% ax2 = subplot(1,2,2); 
%     pcolor(ilon, ilat, phs(1:1:end, 1:1:end)); shading flat; colormap(ax2, flipud(jet)); 
%     h=colorbar; 
%     linkaxes
    


forloopidx = 3:6;
d1 = [];
d2 = []; 
for i=forloopidx
    intdir = cell2mat(intdirs(i)); 
    d1 = [d1; intdir(5:10)];
    d2 = [d2; intdir(12:17)];
    if exist([intfol intdir '/topophase.cor.geo'], 'file')
      % import int 
        filename   = [intfol intdir '/filt_topophase.flat.geo']; 
        fid         = fopen(filename, 'r','native');
        [rmg,count] = fread(fid,[nx*2,ny],'real*4');
        status      = fclose(fid);
        real        = flipud((rmg(1:2:nx*2,1:ny))');
        imag        = flipud((rmg(2:2:nx*2,1:ny))');
        mag         = abs(real+im*imag);
        phs         = angle(real+im*imag);
    end
end
