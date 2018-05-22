d1 = '/data/pmb229/isce/p222f870/mostcombos/int_070713_070828/demLat_N43_N45_Lon_W125_W121.dem'; 
d2 = ''; 

x=importdata([d1 '.vrt']);
l1 = x{1}; 
qf = strfind(l1, '"'); 
nx = str2num(l1(qf(1)+1:qf(2)-1)); 
ny = str2num(l1(qf(3)+1:qf(4)-1)); 

filename2 = [d1]; %OR .FLAT.GEO 
fid          = fopen(filename2,'r','native');
[rmg2,count] = fread(fid,[nx*2,ny],'real*4');