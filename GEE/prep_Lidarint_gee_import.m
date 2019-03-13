% prep_int_gee_import

clear 
close all

%% define files, folders 
pf_fol  = '/data/pmb229/isce/p222f870/'; 
cintfol = 'Lidar_ints/';  
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
    x    = importdata('dl_filt_topophase.unw.geo.vrt');
    l1   = x{1}; 
    qf   = strfind(l1, '"'); 
    nx   = str2num(l1(qf(1)+1:qf(2)-1)); 
    ny   = str2num(l1(qf(3)+1:qf(4)-1)); 

% write a tif file, to get ref frame 
    system('gdal_translate dl_filt_topophase.unw.geo.vrt temptif.tif'); 
    info = geotiffinfo('temptif.tif'); 
cd ..
    

%% separate phs and mag for all ints, write only phs info into geotif
ix       = 0; 
im       = sqrt(-1); 
filename = 'dl_filt_topophase.unw.geo'; 
fid      = fopen(['for_filt_GEE/meta_filt_topophase_unw_geo.csv'], 'wt'); 
fprintf(fid, 'id_no,idx,date1,date2,datenum1,datenum2,baseline\n'); 
for i=1:nints
    bli    = bl(i);
    d1     = datestr(dn(i,1), 'yymmdd');
    d2     = datestr(dn(i,2), 'yymmdd');
    intdir = ['int_' d1 '_' d2]; 
    cd(intdir); 
    
    % import int folder
    fid         = fopen(filename,'r','native');
    [rmg,~]     = fread(fid,[nx,ny*2],'real*4'); 
    mag         = flipud((rmg(1:nx,1:2:ny*2))');
    phs         = flipud((rmg(1:nx,2:2:ny*2))');
%     phs2         = fread(fid,[nx,ny],'real*4'); 
%     phs          = flipud(phs2'); 
    fclose(fid); 

    % find no data pixels, set = -9999
    ndval     = -9999; 
%     midx      = find(mag == 0 & phs == 0);
    midx      = find(phs == 0);
    phs(midx) = ndval;
    
    % write tif file 
    tifname = ['gee_' intdir(5:end) '_filt_topophase_unw_geo.tif']; 
    geotiffwrite(['../for_filt_GEE/' tifname], flipud(phs), info.RefMatrix);
    
    % write csv file 
    fid = fopen(['../for_filt_GEE/meta_filt_topophase_unw_geo.csv'], 'at'); 
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











