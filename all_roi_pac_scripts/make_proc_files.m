% make_proc_files

% ds can be either FIGURE THIS OUT

% function make_proc_files(ds); 

ds = 'baselines/ds_bl-1000m_dl-500m.mat'; 

% variables 
DEMtype='NED'; 
    % MUST put full path to DEM
DEMpath=sprintf('/data/pmb229/roipac/p222f870/DEM_%s/stitched_i2.dem', DEMtype);
filter_strength=0.3; 
range_looks=4; 

ds_exist = exist(ds); 

if ds_exist == 2
    load(ds); 
    ds1 = ds.ds1; 
    ds2 = ds.ds2; 
else
    ds1 = ds(:,1);
    ds2 = ds(:,2);
end

nslc = length(ds1); 

for i = 1:nslc
    s1 = num2str(ds1(i,:)); 
    s2 = num2str(ds2(i,:)); 
    fid = fopen(sprintf('../ints_%s/int_%s_%s_%s.proc', DEMtype, s1, s2, DEMtype), 'wt');
    fprintf(fid, 'SarDir1=../data/%s\n', s1);
    fprintf(fid, 'SarDir2=../data/%s\n', s2);
    fprintf(fid, 'IntDir=int_%s_%s_%s\n', s1, s2, DEMtype);
    fprintf(fid, 'SimDir=int_%s_%s_%s/SIM\n', s1, s2, DEMtype);
    fprintf(fid, 'GeoDir=int_%s_%s_%s/Geo\n', s1, s2, DEMtype);
    fprintf(fid, 'DEM=%s\n', DEMpath); 
    fprintf(fid, 'flattening=topo\n'); 
    fprintf(fid, 'FilterStrength=%g\n', filter_strength); 
    fprintf(fid, 'UnwrappedThreshold = 0.05\n'); 
    fprintf(fid, 'Threshold_mag      = 0.0\n'); 
    fprintf(fid, 'Threshold_ph_grd   = 0.0\n'); 
    fprintf(fid, 'Rlooks_int         = %g\n', range_looks); 
    fprintf(fid, 'Rlooks_unw         = %g\n', range_looks); 
    fprintf(fid, 'Rlooks_geo         = %g\n', range_looks); 
    fprintf(fid, 'pixel_ratio        = 2\n'); 
    fprintf(fid, 'BaslineOrder       = QUAD\n'); 
    fprintf(fid, 'OrbitType          = HDR\n'); 
end










