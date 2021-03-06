%% L1B

cd('/data/pmb229/other/SM/smapL1B/'); 
clear 

d    = dir('SMAP_L1*.h5'); 
d    = {d.name}'; 
info = h5info(cell2mat(d(1)));

% idx of point of interest
p           = [55.105, 19.315];
p    = [55.35, 18.59];

% load data
for i =1:400; %length(d)
    di     = cell2mat(d(i)); 
    hsmh   = hdf5read(di, '/Brightness_Temperature/tb_h_surface_corrected');
    hsmv   = hdf5read(di, '/Brightness_Temperature/tb_v_surface_corrected');
    lati   = hdf5read(di, '/Brightness_Temperature/tb_lat');
    loni   = hdf5read(di, '/Brightness_Temperature/tb_lon');
    
    hsmh(hsmh == -9999) = NaN; 
    hsmv(hsmv == -9999) = NaN; 
    lati(lati == -9999) = NaN; 
    loni(loni == -9999) = NaN; 
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
nn1   = ~isnan(yh); 

%% L1C

cd('/data/pmb229/other/SM/smapL1C/'); 

d    = dir('SMAP_L1*.h5'); 
d    = {d.name}'; 
info = h5info(cell2mat(d(1)));

% idx of point of interest
p           = [55.105, 19.315];
p    = [55.35, 18.59];

% load data
for i =1:400 %length(d)
    di     = cell2mat(d(i)); 
    hsmh   = hdf5read(di, '/Global_Projection/cell_tb_h_surface_corrected_aft');
    hsmv   = hdf5read(di, '/Global_Projection/cell_tb_v_surface_corrected_aft');
    errh   = hdf5read(di, '/Global_Projection/cell_tb_error_h_aft'); 
    errv   = hdf5read(di, '/Global_Projection/cell_tb_error_v_aft'); 
    lati   = hdf5read(di, '/Global_Projection/cell_lat_centroid_aft');
    loni   = hdf5read(di, '/Global_Projection/cell_lon_centroid_aft');
    
    hsmh(hsmh == -9999) = NaN; 
    hsmv(hsmv == -9999) = NaN; 
    lati(lati == -9999) = NaN; 
    loni(loni == -9999) = NaN; 
    dt2(i,:)            = di(21:28);
    ad(i,:)             = di(19); 
    
    ll = [loni(:) lati(:)]; 
    [idx, dist] = dsearchn(ll,p); 
    if dist < 0.35
        [ii,jj] = ind2sub(size(lati), idx); 
        yh2(i) = hsmh(ii,jj); 
        yv2(i) = hsmv(ii,jj); 
        eh2(i) = errh(ii,jj); 
        ev2(i) = errv(ii,jj); 
    else
        yh2(i) = NaN; 
        yv2(i) = NaN; 
        eh2(i) = NaN; 
        ev2(i) = NaN; 
    end
end
dt2                = datenum(dt2, 'yyyymmdd');
nn2 = ~isnan(yh2); 

%% L2

cd('/data/pmb229/other/SM/smapL2/'); 


d    = dir('SMAP_L2*.h5'); 
d    = {d.name}'; 
info = h5info(cell2mat(d(1)));

% idx of point of interest
%p   = [55.105, 19.315];
p    = [54.4 17.9]; 
p    = [55.35, 18.59];

% load data
for i =1:460 %length(d)
    di     = cell2mat(d(i)); 
    hsm    = hdf5read(di, '/Soil_Moisture_Retrieval_Data/soil_moisture_option1'); % option1 = H-pol (best), 2=v-pol, 3=dual pole
    lati   = hdf5read(di, '/Soil_Moisture_Retrieval_Data/latitude');
    loni   = hdf5read(di, '/Soil_Moisture_Retrieval_Data/longitude');
     dt3(i,:)   = di(22:29);
    
    ll = [loni lati]; 
    [idx, dist] = dsearchn(ll,p); 
    if dist < 0.35
        [ii,jj] = ind2sub(size(lati), idx); 
        y3(i) = hsm(ii,jj); 
    else
        y3(i) = NaN; 
    end
end
dt3                = datenum(dt3, 'yyyymmdd');
nn3 = ~isnan(y3); 

%% plot point of interest
close all; 

fig = figure('units', 'normalized', 'outerposition', [.1 .7 .7 .5]);
left_color = [.0 .0 1];
right_color = [1 1 1]*0;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
hold on;  box on; 

%errorbar(dt, yh, ones(1, length(yh))*1.3, 'r.', 'markersize', 10); 
%errorbar(dt, yv, ones(1, length(yh))*1.3, 'r.', 'markersize', 10); 
% errorbar(dt2, yh2, eh2, 'r.', 'markersize', 10); 
% errorbar(dt2, yv2, ev2, 'b.', 'markersize', 10); 

yyaxis left
plot(dt(nn1), yh(nn1), 'bo--', 'markersize', 6); 
plot(dt(nn1), yv(nn1)-30, 'b.--', 'markersize', 15); 
plot(dt2(nn2), yh2(nn2), 'bo-', 'markersize', 6); 
plot(dt2(nn2), yv2(nn2)-30, 'b.-', 'markersize', 15); 
ylabel('Brightness Temp'); 

yyaxis right
plot(dt3(nn3), 1-y3(nn3), 'k.-', 'markersize', 15); 
ylabel('Soil Moisture'); 

 
xlabel('Date'); 
xlim([datenum('2018-05-01') datenum('2018-06-15')]);
datetick('x', 'keeplimits');

legend('L1B (h)', 'L1B (v)', 'L1C (h)', 'L1C (v)', 'L2', 'location', 'southwest');
% legend('L1C (h)', 'L1C (v)', 'L2', 'location', 'southwest');
dcmObj = datacursormode;
set(dcmObj, 'UpdateFcn', @datacursorprecision);














