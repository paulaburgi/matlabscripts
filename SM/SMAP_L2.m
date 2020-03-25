cd('/data/pmb229/other/SM/smapL2/'); 
clear 

d    = dir('SMAP_L2*.h5'); 
d    = {d.name}'; 
info = h5info(cell2mat(d(1)));

% idx of point of interest
%p   = [55.105, 19.315];
p    = [54.4 17.9]; 
p    = [55.35, 18.59];

% load data
n = length(d);
sm  = cell(n,1); 
lat = cell(n,1); 
lon = cell(n,1); 
for i =1:length(d)
    di     = cell2mat(d(i)); 
    hsm    = hdf5read(di, '/Soil_Moisture_Retrieval_Data/soil_moisture_option1'); % option1 = H-pol (best), 2=v-pol, 3=dual pole
    lati   = hdf5read(di, '/Soil_Moisture_Retrieval_Data/latitude');
    loni   = hdf5read(di, '/Soil_Moisture_Retrieval_Data/longitude');
    
    sm(i)     = {hsm(hsm == -9999)}; 
    lat(i)    = {lati(lati == -9999)}; 
    lon(i)    = {loni(loni == -9999)}; 
    dt(i,:)   = di(22:29);
    
    ll = [loni lati]; 
    [idx, dist] = dsearchn(ll,p); 
    if dist < 0.35
        [ii,jj] = ind2sub(size(lati), idx); 
        y(i) = hsm(ii,jj); 
    else
        y(i) = NaN; 
    end
end
dt                = datenum(dt, 'yyyymmdd');


%% plot point of interest
%close all; 

yn = ~isnan(y); 


%figure('units', 'normalized', 'outerposition', [.1 .7 .9 .5]);
hold on; grid minor; box on; 
plot(dt(yn), y(yn), '.-', 'markersize', 10); 
datetick; 
xlabel('date'); 
ylabel('soil moisture'); 
dcmObj = datacursormode;
set(dcmObj, 'UpdateFcn', @datacursorprecision);















