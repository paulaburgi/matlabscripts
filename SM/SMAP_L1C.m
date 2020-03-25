cd('/data/pmb229/other/SM/smapL1C/'); 
clear 

d    = dir('SMAP_L1*.h5'); 
d    = {d.name}'; 
info = h5info(cell2mat(d(1)));

% idx of point of interest
p           = [55.105, 19.315];

% load data
for i =1:length(d)
    di     = cell2mat(d(i)); 
    hsmh   = hdf5read(di, '/Global_Projection/cell_tb_h_surface_corrected_aft');
    hsmv   = hdf5read(di, '/Global_Projection/cell_tb_h_surface_corrected_fore');
    lati   = hdf5read(di, '/Global_Projection/cell_lat_centroid_aft');
    loni   = hdf5read(di, '/Global_Projection/cell_lon_centroid_aft');
    
    hsmh(hsmh == -9999) = NaN; 
    hsmv(hsmv == -9999) = NaN; 
    lati(lati == -9999) = NaN; 
    loni(loni == -9999) = NaN; 
    h(i)                = {hsmh}; 
    v(i)                = {hsmh}; 
    lat(i)              = {lati}; 
    lon(i)              = {loni}; 
    dt(i,:)             = di(21:28);
    ad(i,:)             = di(19); 
    
    ll = [loni(:) lati(:)]; 
    [idx, dist] = dsearchn(ll,p); 
    if dist < 0.35
        [ii,jj] = ind2sub(size(lati), idx); 
        yh(i) = hsmh(ii,jj); 
        yv(i) = hsmv(ii,jj); 
    else
        yh(i) = NaN; 
        yv(i) = NaN; 
    end
end
dt                = datenum(dt, 'yyyymmdd');


%% plot point of interest
%close all; 

figure; hold on; grid minor; box on; 
for i = 1:length(d)
    if strcmp(ad(i), 'A')
        plot(dt(i), yh(i), 'ro', 'markersize', 10); 
        plot(dt(i), yv(i), 'bo', 'markersize', 10); 
    else
        plot(dt(i), yh(i), 'rx', 'markersize', 10); 
        plot(dt(i), yv(i), 'bx', 'markersize', 10); 
    end
end

datetick; 
xlabel('date'); 
ylabel('brightness temp'); 
legend('h', 'v');















