% prep_int_gee_import

clear 
close all

%% define files, folders 
% sumatra
% pf_fol  = '/data/pmb229/isce/p446f7190_sumatra/'; 
% cintfol = 'ints_SRTM/';  pol = 'HH'; 
% datafol = [pf_fol 'data/']; 
% intfol  = [pf_fol cintfol];

% cascadia
pf_fol  = '/data/pmb229/isce/p222f870/'; 
cintfol = 'mostcombos/';  pol = 'HH'; 
datafol = [pf_fol 'data/']; 
intfol  = [pf_fol cintfol];

cd(intfol); 
intdirs = dir; 
intdirs = {intdirs.name}; 
load([datafol 'analysis/meancor_bl_dates_area2_' pol '.mat']); 
dc      = meancor_bl_dates.dateCombos; 

%% get int size, write temp geotiff
cd(cell2mat(intdirs(3))); 
% get nx, ny
    x  = importdata('filt_topophase.flat.geo.vrt');
    l1 = x{1}; 
    qf = strfind(l1, '"'); 
    nx = str2num(l1(qf(1)+1:qf(2)-1)); 
    ny = str2num(l1(qf(3)+1:qf(4)-1)); 

% write a tif file, to get ref frame 
    system('gdal_translate filt_topophase.flat.geo.vrt temptif.tif'); 
    info = geotiffinfo('temptif.tif'); 
cd ..
    

%% separate phs and mag for all ints, write only phs info into geotif
ix=-1; 
for i=3:length(intdirs); 
    intdir = cell2mat(intdirs(i)); 
    cd(intdir);
    d1 = [intdir(5:10)];
    d2 = [intdir(12:17)];
    dn = [datenum(d1, 'yymmdd') datenum(d2, 'yymmdd')]; 
    ix = ix+1; 
    
    deq1 = eq(dc, dn); 
    didx  = find(deq1(:,1) == 1 & deq1(:,2) == 1); 
    
    % import int folder
    im = sqrt(-1); 
    filename   = 'filt_topophase.flat.geo'; 
    fid         = fopen(filename, 'r','native');
    [rmg,count] = fread(fid,[nx*2,ny],'real*4');
    status      = fclose(fid);
    real        = flipud((rmg(1:2:nx*2,1:ny))');
    imag        = flipud((rmg(2:2:nx*2,1:ny))');
    mag         = abs(real+im*imag);
    phs         = angle(real+im*imag);
    
    % get baseline info
    bl = meancor_bl_dates.bl(didx); 
    
    % find no data pixels, set = -9999
    ndval    = -9999; 
    midx      = find(mag == 0 & phs == 0);
    phs(midx) = ndval;
    
    % write tif file 
    tifname = ['gee_' intdir(5:end) '_filt_topophase_flat_geo.tif']; 
    geotiffwrite(tifname, flipud(phs), info.RefMatrix);
    
    % write csv file 
    fid = fopen(['meta_' filename '.csv'], 'wt'); 
    fprintf(fid, 'id_no,idx,date1,date2,datenum1,datenum2,baseline\n'); 
    fprintf(fid, [tifname(1:end-4) ',']); 
    fprintf(fid, [num2str(ix) ',']); 
    fprintf(fid, ['d' d1 ',']); 
    fprintf(fid, ['d' d2 ',']); 
    fprintf(fid, [num2str(dn(1)) ',']); 
    fprintf(fid, [num2str(dn(2)) ',']); 
    fprintf(fid, num2str(round(bl))); 
    fclose(fid); 
    
    cd ..
end











