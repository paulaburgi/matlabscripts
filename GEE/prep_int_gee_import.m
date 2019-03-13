% prep_int_gee_import

clear 
close all

%% define files, folders 
% sumatra
pf_fol  = '/data/pmb229/isce/p446f7190_sumatra/'; 
cintfol = 'ints/';  
datafol = [pf_fol 'data/']; 
intfol  = [pf_fol cintfol];

% cascadia
% pf_fol  = '/data/pmb229/isce/p222f870/'; 
% cintfol = 'mostcombos/';  pol = 'HH'; 
% datafol = [pf_fol 'data/']; 
% intfol  = [pf_fol cintfol];

cd(intfol); 
intdirs = dir('int_*'); 
intdirs = {intdirs.name}; 


%% get int size, write temp geotiff
cd(cell2mat(intdirs(1))); 
% get nx, ny
    x  = importdata('filt_topophase.unw.geo.vrt');
    l1 = x{1}; 
    qf = strfind(l1, '"'); 
    nx = str2num(l1(qf(1)+1:qf(2)-1)); 
    ny = str2num(l1(qf(3)+1:qf(4)-1)); 

% write a tif file, to get ref frame 
    system('gdal_translate filt_topophase.unw.geo.vrt tempunwtif.tif'); 
    info = geotiffinfo('tempunwtif.tif'); 
cd ..

geefold = '/data/pmb229/isce/p446f7190_sumatra/ints/for_gee/'; 
fid = fopen([geefold 'meta_int_all.csv'], 'wt'); 
fprintf(fid, 'id_no,idx,date1,date2,datenum1,datenum2,baseline\n'); 

%% separate phs and mag for all ints, write only phs info into geotif
ix=-1; 
for i=1:length(intdirs); 
    intdir = cell2mat(intdirs(i)); 
    cd(intdir);
    d1 = [intdir(5:10)];
    d2 = [intdir(12:17)];
    dn = [datenum(d1, 'yymmdd') datenum(d2, 'yymmdd')]; 
    ix = ix+1; 
    
    % temp tiff 
    % get nx, ny
    x  = importdata('filt_topophase.unw.geo.vrt');
    l1 = x{1}; 
    qf = strfind(l1, '"'); 
    nx = str2num(l1(qf(1)+1:qf(2)-1)); 
    ny = str2num(l1(qf(3)+1:qf(4)-1)); 

% write a tif file, to get ref frame 
    system('gdal_translate filt_topophase.unw.geo.vrt tempunwtif.tif'); 
    info = geotiffinfo('tempunwtif.tif'); 
    
    % baseline
    [~, t] = system('more isce.log | grep perp_baseline_top');
    [~, b] = system('more isce.log | grep perp_baseline_bottom');
    bl = round(mean([str2num(t(strfind(t, '=')+1:end)) str2num(b(strfind(b, '=')+1:end))])); 
    
    % import int folder
    filename    = 'filt_topophase.unw.geo'; 
    fid         = fopen(filename,'r','native');
    [rmg,~]     = fread(fid,[nx,ny*2],'real*4'); 
                  fclose(fid); 
    mag         = flipud((rmg(1:nx,1:2:ny*2))');
    unw         = flipud((rmg(1:nx,2:2:ny*2))');
    
    % find no data pixels, set = -9999
    ndval    = -9999; 
    midx      = find(mag == 0 & unw == 0);
    unw(midx) = ndval;
    
    % write tif file 
    tifname = ['gee_' intdir(5:end) '_unw.tif']; 
    geotiffwrite([geefold tifname], flipud(unw), info.RefMatrix);
    
    
    % write csv file 
    fid = fopen([geefold 'meta_int_all.csv'], 'a'); 
    fprintf(fid, [tifname(1:end-4) ',']); 
    fprintf(fid, [num2str(ix) ',']); 
    fprintf(fid, ['d' d1 ',']); 
    fprintf(fid, ['d' d2 ',']); 
    fprintf(fid, [num2str(dn(1)) ',']); 
    fprintf(fid, [num2str(dn(2)) ',']); 
    fprintf(fid, [num2str(round(bl)) '\n']); 
    fclose(fid); 
    
    cd ..
end











