% make .proc files for ALOS data

% THIS IS USELESS. MAKE SAME SCRIPT, EXCEPT INSTEAD OF SEQUENTIAL, MAKE IT
% BE A MIN BASELINE REQUIREMENT

clear
s = dir; 

% variables 
DEMtype='NED'; 
    % MUST put full path to DEM
DEMpath='/data/pmb229/isce/p447f7170_sumatra/dem/p447f7170.dem';
filter_strength=0.3; 
range_looks=4; 


% find file/folders with digits in file/folder name
f=[]; 
for i=1:length(s)
    TF = isstrprop(s(i).name,'digit'); 
    f = [f; sum(TF)];
end

% find indices of files slc folders 
slcidx = find(f==6);
nslc   = length(slcidx);
% test
% concatenate all slc folder names
fids=[]; 
for i=1:nslc
    sn = s(slcidx(i)).name; 
    fids=[fids; sn];
end

s1 = fids(1,:); 
% make .proc files 

bld = 'baselines';
ble = exist(bld); 

if ble == 7 
    cd(bld);
else
    mkdir(bld);
    cd(bld); 
end


for i = 2:nslc
    s2 = fids(i,:); 
    fid = fopen(sprintf('int_%s_%s_%s.proc', s1, s2, DEMtype), 'wt');
    fprintf(fid, 'SarDir1=../%s\n', s1);
    fprintf(fid, 'SarDir2=../%s\n', s2);
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


