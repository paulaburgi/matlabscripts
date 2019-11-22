files=dir('*/*sm_rootzone*tif');
dfile=files(1:2:end);
efile=files(2:2:end);

n=length(dfile);

for i=1:n
    name=dfile(i).name;
    dn(i)=datenum(name(16:26),'YYYYmmDDTHH');
    dat=geotiffread([dfile(i).folder '/' name]);
    alldat(i,:,:)=dat;
     dat=geotiffread([efile(i).folder '/' efile(i).name]);
   allerr(i,:,:)=dat;
end






d    = dir('SMAP_L3*.h5'); 
d    = {d.name}'; 
info = h5info(cell2mat(d(1))); 

p      = [55.105, 19.315];
x      = [];
y      = [];
x_pm   = [];
y_pm   = [];
for i =1:length(d)
    di     = cell2mat(d(i)); 
    dt     = datenum(di(14:21), 'yyyymmdd');
    hlonc  = hdf5read(di, '/Soil_Moisture_Retrieval_Data_AM/longitude');
    hlatc  = hdf5read(di, '/Soil_Moisture_Retrieval_Data_AM/latitude');
    hsm    = hdf5read(di, '/Soil_Moisture_Retrieval_Data_AM/soil_moisture');
    hloncp = hdf5read(di, '/Soil_Moisture_Retrieval_Data_PM/longitude_pm');
    hlatcp = hdf5read(di, '/Soil_Moisture_Retrieval_Data_PM/latitude_pm');
    hsmp   = hdf5read(di, '/Soil_Moisture_Retrieval_Data_PM/soil_moisture_pm');
    
    % am
    ll          = [hlonc(:) hlatc(:)]; 
    [idx, dist] = dsearchn(ll, p); 
    if dist < 0.35
        si = hsm(:); 
        y  = [y; si(idx)]; 
    else
        y  = [y; NaN]; 
    end
    x  = [x; double(dt) ll(idx,:)];
    % pm
    llp          = [hloncp(:) hlatcp(:)]; 
    [idxp, distp] = dsearchn(llp, p); 
    if distp < 0.35
        sip = hsmp(:); 
        y_pm  = [y_pm; sip(idxp)]; 
    else
        y_pm  = [y_pm; NaN]; 
    end
    x_pm  = [x_pm; double(dt) ll(idxp,:)];

end