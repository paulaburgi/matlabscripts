% prep_int_gee_import

clear 
close all

%% define files, folders 
pf_fol  = '/data/pmb229/isce/p222f870/'; 
cintfol = 'NED_ints/';  
datafol = [pf_fol 'data/']; 
intfol  = [pf_fol cintfol];

cd(intfol); 

load('log/useints.mat'); 
dn    = useints.dn; 
bl    = useints.bl; 
nints = length(dn); 


%% get int size, write temp geotiff
intdir1  = ['int_' datestr(dn(1,1), 'yymmdd') '_' datestr(dn(1,2), 'yymmdd')]; 
cd(intdir1);
% get nx, ny
    x    = importdata('dl_filt_topophase.flat.geo.vrt');
    l1   = x{1}; 
    qf   = strfind(l1, '"'); 
    nx   = str2num(l1(qf(1)+1:qf(2)-1)); 
    ny   = str2num(l1(qf(3)+1:qf(4)-1)); 

% write a tif file, to get ref frame 
    system('gdal_translate dl_filt_topophase.flat.geo.vrt temptif.tif'); 
    info = geotiffinfo('temptif.tif'); 
cd ..
    

%% separate phs and mag for all ints, write only phs info into geotif
ix       = 0; 
im       = sqrt(-1); 
filename = 'dl_filt_topophase.flat.geo'; 
fid      = fopen(['for_GEE2/meta_filt_topophase_flat_geo.csv'], 'wt'); 
fprintf(fid, 'id_no,idx,date1,date2,datenum1,datenum2,baseline\n'); 
for i=1:nints
    bli    = bl(i);
    d1     = datestr(dn(i,1), 'yymmdd');
    d2     = datestr(dn(i,2), 'yymmdd');
    intdir = ['int_' d1 '_' d2]; 
    cd(intdir); 
    
    % import int folder
    fid         = fopen(filename, 'r','native');
    [rmg,count] = fread(fid,[nx*2,ny],'real*4');
    status      = fclose(fid);
    real        = flipud((rmg(1:2:nx*2,1:ny))');
    imag        = flipud((rmg(2:2:nx*2,1:ny))');
    mag         = abs(real+im*imag);
    phs         = angle(real+im*imag);
    
    % find no data pixels, set = -9999
    ndval    = -9999; 
    midx      = find(mag == 0 & phs == 0);
    phs(midx) = ndval;
    
    % write tif file 
    tifname = ['gee_' intdir(5:end) '_filt_topophase_flat_geo.tif']; 
    geotiffwrite(['../for_GEE2/' tifname], flipud(phs), info.RefMatrix);
    
    % write csv file 
    fid = fopen(['../for_GEE2/meta_filt_topophase_flat_geo.csv'], 'at'); 
    fprintf(fid, [tifname(1:end-4) ',']); 
    fprintf(fid, [num2str(ix) ',']); 
    fprintf(fid, ['d' d1 ',']); 
    fprintf(fid, ['d' d2 ',']); 
    fprintf(fid, [num2str(dn(i,1)) ',']); 
    fprintf(fid, [num2str(dn(i,2)) ',']); 
    fprintf(fid, [num2str(round(bli)) '\n']); 
    fclose(fid); 
    
    cd ..
    ix = ix+1; 
end











