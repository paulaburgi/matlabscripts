% prep_cor_gee_import

clear 
close all

%% define files, folders 
% cascadia
% pf_fol  = '/data/pmb229/isce/p222f870/'; 
% cintfol = '/data/pmb229/isce/p222f870/Lidar_ints/';  pol = 'HH'; 
% datafol = [pf_fol 'data/']; 
% intfol  = [cintfol];

% sumatra
pf_fol  = '/data/pmb229/isce/p446f7190_sumatra/'; 
cintfol = 'ints/';  
datafol = [pf_fol 'data/']; 
intfol  = [pf_fol cintfol];

cd(intfol); 

% for srtm ints
intdirs = dir('int_*'); 
intdirs = {intdirs.name}; 
% load([datafol 'analysis/meancor_bl_dates_area2_' pol '.mat']); 
% dc      = meancor_bl_dates.dateCombos; 
nints = length(intdirs);

% for lidar ints
% load('log/useints.mat'); 
% dn    = useints.dn; 
% bl    = useints.bl; 
% nints = length(dn); 

%% get int size, write temp geotiff
% for srtm ints
    % cd(cell2mat(intdirs(1))); 
% for lidar ints
%     intdir1  = ['int_' datestr(dn(1,1), 'yymmdd') '_' datestr(dn(1,2), 'yymmdd')]; 
%     cd(intdir1);
% get nx, ny
%     x  = importdata('topophase.cor.geo.vrt');
%     l1 = x{1}; 
%     qf = strfind(l1, '"'); 
%     nx = str2num(l1(qf(1)+1:qf(2)-1)); 
%     ny = str2num(l1(qf(3)+1:qf(4)-1)); 
% 
% % write a tif file, to get ref frame 
%     system('gdal_translate topophase.cor.geo.vrt tempcortif.tif'); 
%     info = geotiffinfo('temptif.tif'); 
% cd ..
    

%% separate cor and mag for all ints, write only cor info into geotif
ix=-1; 
fid      = fopen(['for_cor_GEE/meta_cor_geo.csv'], 'wt'); 
fprintf(fid, 'id_no,idx,date1,date2,datenum1,datenum2,baseline\n'); 
dn = []; 
for i=1:nints
    % for srtm
       intdir = cell2mat(intdirs(i)); 
       cd(intdir);
       d1 = [intdir(5:10)];
        d2 = [intdir(12:17)];
        dn = [dn; [datenum(d1, 'yymmdd') datenum(d2, 'yymmdd')]];
%         deq1 = eq(dc, dn); 
%         didx  = find(deq1(:,1) == 1 & deq1(:,2) == 1); 

% for lidar
%         bli    = bl(i);
%         d1     = datestr(dn(i,1), 'yymmdd');
%         d2     = datestr(dn(i,2), 'yymmdd');
%         intdir = ['int_' d1 '_' d2]; 
%         cd(intdir); 

    ix = ix+1; 
    
    % get nx, ny
    x  = importdata('topophase.cor.geo.vrt');
    l1 = x{1}; 
    qf = strfind(l1, '"'); 
    nx = str2num(l1(qf(1)+1:qf(2)-1)); 
    ny = str2num(l1(qf(3)+1:qf(4)-1)); 

% write a tif file, to get ref frame 
    system('gdal_translate topophase.cor.geo.vrt tempcortif.tif'); 
    info = geotiffinfo('tempcortif.tif'); 
    
    % import int folder
    filename    = 'topophase.cor.geo'; 
    fid         = fopen(filename,'r','native');
    [rmg,~]     = fread(fid,[nx,ny*2],'real*4'); 
                  fclose(fid); 
    mag         = flipud((rmg(1:nx,1:2:ny*2))');
    cor         = flipud((rmg(1:nx,2:2:ny*2))');
    
    % get baseline info
%     bl = meancor_bl_dates.bl(didx); 
    % baseline
    [~, t] = system('more isce.log | grep perp_baseline_top');
    [~, b] = system('more isce.log | grep perp_baseline_bottom');
    bli = round(mean([str2num(t(strfind(t, '=')+1:end)) str2num(b(strfind(b, '=')+1:end))])); 

    
    % find no data pixels, set = -9999
    ndval    = -9999; 
    midx      = find(mag == 0 & cor == 0);
    cor(midx) = ndval;
    
    % write tif file 
    tifname = ['gee_' intdir(5:end) '_topophase_cor_geo.tif']; 
    geotiffwrite(tifname, flipud(cor), info.RefMatrix);
    system(['mv ' tifname ' ../for_cor_GEE']); 
    
    % write csv file 
    fid = fopen(['../for_cor_GEE/meta_cor_geo.csv'], 'at'); 
    fprintf(fid, [tifname(1:end-4) ',']); 
    fprintf(fid, [num2str(ix) ',']); 
    fprintf(fid, ['d' d1 ',']); 
    fprintf(fid, ['d' d2 ',']); 
    fprintf(fid, [num2str(dn(i,1)) ',']); 
    fprintf(fid, [num2str(dn(i,2)) ',']); 
    fprintf(fid, [num2str(round(bli)) '\n']); 
    fclose(fid); 
    
    cd ..
end











