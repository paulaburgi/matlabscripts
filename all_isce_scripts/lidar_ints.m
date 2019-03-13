cd('/data/pmb229/isce/p222f870/Lidar_ints/'); 

load('log/useints.mat'); 
dn = useints.dn; 
bl = useints.bl; 
nints = length(dn); 


% link IMG, LED, dem
% cp over int.xml files
% replace dem with lidar dem 
% run all 

for i = 1:nints
    d1 = datestr(dn(i,1), 'yymmdd'); 
    d2 = datestr(dn(i,2), 'yymmdd'); 
    intdir = ['int_' d1 '_' d2]; 
    
    % make directory
    if ~exist(intdir, 'file')
        system(['mkdir ' intdir]); 
    end
    cd(intdir); 
    
    % link IMG, LED
    dd1 = ['20' d1]; 
    dd2 = ['20' d2]; 
    system(['ln -sf /data/pmb229/isce/p222f870/data/' dd1 '/IMG-HH*']); 
    system(['ln -sf /data/pmb229/isce/p222f870/data/' dd2 '/IMG-HH*']); 
    system(['ln -sf /data/pmb229/isce/p222f870/data/' dd1 '/LED*']); 
    system(['ln -sf /data/pmb229/isce/p222f870/data/' dd2 '/LED*']); 
    
    % link dem
    system('ln -sf /data/pmb229/isce/p222f870/DEMs/OregonLidar/OregonLidar_median.dem'); 
    system('ln -sf /data/pmb229/isce/p222f870/DEMs/OregonLidar/OregonLidar_median.dem.vrt'); 
    system('ln -sf /data/pmb229/isce/p222f870/DEMs/OregonLidar/OregonLidar_median.dem.xml'); 
    
    % cp xml files from NED ints
    neddir = ['/data/pmb229/isce/p222f870/NED_ints/' intdir '/' intdir '.xml'];
    system(['cp ' neddir ' .']); 
    
    % replace ned dem with lidar dem
    old = 'stitched_i2.dem.xml'; 
    new = 'OregonLidar_median.dem.xml'; 
    system(['sed -i ''s/' old '/' new '/g'' ' intdir '.xml']);
    
    % set so it unwraps
    old = 'False'; 
    new = 'True'; 
    system(['sed -i ''s/' old '/' new '/g'' ' intdir '.xml']);
    old = 'grass'; 
    new = 'snaphu'; 
    system(['sed -i ''s/' old '/' new '/g'' ' intdir '.xml']);
    
    % geocode unwrapped product
    old = '"filt_topophase.flat", "topophase.flat"'; 
    new = '"filt_topophase.unw"'; 
    system(['sed -i ''s/' old '/' new '/g'' ' intdir '.xml']);
    
    % geocode unwrapped product
    old = '"filt_topophase.unw"'; 
    new = '"topophase.cor"'; 
    system(['sed -i ''s/' old '/' new '/g'' ' intdir '.xml']);
    
    % change filtering
    old = '0.3'; 
    new = '0.5'; 
    system(['sed -i ''s/' old '/' new '/g'' ' intdir '.xml']);
    
    cd ..
end

