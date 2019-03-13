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
cd('geocodedSLCs/'); 
slcdirs = dir('20*geo'); 
slcdirs = {slcdirs.name}'; 

%% get slc size, write temp geotiff
temp = cell2mat(slcdirs(1)); 
% get nx, ny
    x  = importdata([temp '.vrt']);
    l1 = x{1}; 
    qf = strfind(l1, '"'); 
    nx = str2num(l1(qf(1)+1:qf(2)-1)); 
    ny = str2num(l1(qf(3)+1:qf(4)-1)); 

% write a tif file, to get ref frame 
    system(['gdal_translate ' temp '.vrt temptif.tif']); 
    info = geotiffinfo('temptif.tif'); 
    
geefold = './for_gee/'; 
fid = fopen([geefold 'meta_slc_all.csv'], 'wt'); 
fprintf(fid, 'id_no,idx,date\n'); 

%% separate phs and mag for all ints, write only phs info into geotif
ix=-1; 
for i=1:length(slcdirs); 
    slcdir = cell2mat(slcdirs(i)); 
    d = [slcdir(1:8)];
    ix = ix+1; 
  
    % import int folder
    im = sqrt(-1); 
    filename   = slcdir; 
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
    mag(midx) = ndval;
    
    % write tif file
    tifname = ['gee_slcmag_' slcdir(1:8) '.tif']; 
    geotiffwrite([geefold tifname], flipud(mag), info.RefMatrix);
    
    % write csv file 
    fid = fopen([geefold 'meta_slc_all.csv'], 'a'); 
    fprintf(fid, [tifname(1:end-4) ',']); 
    fprintf(fid, [num2str(ix) ',']); 
    fprintf(fid, [d '\n']); 
    fclose(fid); 
    
end











