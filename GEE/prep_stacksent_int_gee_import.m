% prep_int_gee_import

clear 
close all

%% define files, folders 
% cascadia
pf_fol  = '/data/pmb229/isce/p222f870/Sentinel_p13f445/'; 
cintfol = 'lidar_ints/merged/interferograms/';  
datafol = [pf_fol 'data/']; 
intfol  = [pf_fol cintfol];

cd(intfol); 
intdirs = dir('2015*'); 
intdirs = {intdirs.name}; 


%% get int size, write temp geotiff
cd(cell2mat(intdirs(1))); 
% get nx, ny
    x  = importdata('filt_fine.unw.geo.vrt');
    l1 = x{1}; 
    qf = strfind(l1, '"'); 
    nx = str2num(l1(qf(1)+1:qf(2)-1)); 
    ny = str2num(l1(qf(3)+1:qf(4)-1)); 

% write a tif file, to get ref frame 
    system('gdal_translate filt_fine.unw.geo.vrt tempunwtif.tif'); 
    info = geotiffinfo('tempunwtif.tif'); 
cd ..

geefold = '/data/pmb229/isce/p222f870/Sentinel_p13f445/lidar_ints/for_gee/'; 
fid = fopen([geefold 'meta_int_all.csv'], 'wt'); 
fprintf(fid, 'id_no,idx,date1,date2,datenum1,datenum2,baseline\n'); 

bdir  = '/data/pmb229/isce/p222f870/Sentinel_p13f445/lidar_ints/baselines/'; 
bldir = dir([bdir '2015*']); 
bldir = {bldir.name}';
bla   = 0; 
dt1   = cell2mat(bldir(1));
dta   = {dt1(1:8)}; 
for i = 1:length(bldir)
    di      = cell2mat(bldir(i)); 
    [~, t]  = system(['more ' bdir di '/*.txt | grep Bperp']); 
    ci      = strfind(t, ':'); 
    bl      = round(mean([str2num(t(ci(1)+1:ci(1)+5)) str2num(t(ci(2)+1:ci(2)+5))])); 
    bla     = [bla; bl]; 
    dta{i+1,1}     = di(10:17);
end

%% separate phs and mag for all ints, write only phs info into geotif
ix=-1; 
for i=1:length(intdirs); 
    intdir = cell2mat(intdirs(i)); 
    cd(intdir);
    d1 = [intdir(1:8)];
    d2 = [intdir(10:17)];
    dn = [datenum(d1, 'yyyymmdd') datenum(d2, 'yyyymmdd')]; 
    ix = ix+1; 
    
    % baseline
    bidx1 = bla(find(strcmp(dta, d1)));
    bidx2 = bla(find(strcmp(dta, d2)));
    bl    = bidx2-bidx1; 
    
    
    % import int folder
    filename    = 'filt_fine.unw.geo'; 
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
    tifname = ['gee_' intdir '_unw.tif']; 
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











