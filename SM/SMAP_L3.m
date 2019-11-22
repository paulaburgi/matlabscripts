cd('/data/pmb229/other/smapL3/'); 
clear 

d    = dir('SMAP_L3*.h5'); 
d    = {d.name}'; 
info = h5info(cell2mat(d(1)));

% load data
for i =1:length(d)
    di     = cell2mat(d(i)); 
    hsm    = hdf5read(di, '/Soil_Moisture_Retrieval_Data_AM/soil_moisture');
    hsmp   = hdf5read(di, '/Soil_Moisture_Retrieval_Data_PM/soil_moisture_pm');
    
    sm_am(i,:,:) = hsm; 
    sm_pm(i,:,:) = hsmp; 
    dt(i,:)      = di(14:21);
end
sm_am(find(sm_am == -9999)) = NaN; 
sm_pm(find(sm_pm == -9999)) = NaN; 
dt                          = datenum(dt, 'yyyymmdd');

% idx of point of interest
p           = [55.105, 19.315];
hlonc       = hdf5read(cell2mat(d(1)), '/Soil_Moisture_Retrieval_Data_AM/longitude');
hlatc       = hdf5read(cell2mat(d(1)), '/Soil_Moisture_Retrieval_Data_AM/latitude');
ll          = [hlonc(:) hlatc(:)]; 
[idx, dist] = dsearchn(ll, p); % check dist is <0.35

%% plot point of interest
close all; 
[i,j] = ind2sub(size(hlonc), idx); 
y_am  = sm_am(:,i,j); 
y_pm  = sm_pm(:,i,j); 

figure; hold on; grid minor; box on; 
plot(dt, y_am, 'r.', 'markersize', 10); 
plot(dt, y_pm, 'b.', 'markersize', 10); 
ylim([min([y_am; y_pm])-0.01  max([y_am; y_pm])+0.01]); 
xlim([min(dt)-10              max(dt)+10]);
datetick; 
xlabel('date'); 
ylabel('soil moisture'); 
legend('am', 'pm');




















