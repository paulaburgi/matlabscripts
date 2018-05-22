% concat_gee_csvs_tifs

clear 
close all

%% define files, folders 
% sumatra
pf_fol  = '/data/pmb229/isce/p446f7190_sumatra/'; 
cintfol = 'ints_SRTM/';  pol = 'HH'; 

% cascadia
pf_fol  = '/data/pmb229/isce/p222f870/'; 
cintfol = 'mostcombos/';  pol = 'HH'; 

datafol = [pf_fol 'data/']; 
intfol  = [pf_fol cintfol];
geefol = [datafol 'analysis/geotiff_gee/']; 


cd(intfol);  
csv_all = [geefol 'meta_all.csv'];

% only concatenate csvs for ints that have good correlation
load([datafol 'analysis/meancor_bl_dates_area2_' pol '.mat']); 
    datesall  = meancor_bl_dates.dateCombos;     
    idx       = meancor_bl_dates.good_cor_idx;   
    dates     = datesall(idx,:); 
    d1_all    = datestr(dates(:,1), 'yymmdd'); 
    d2_all    = datestr(dates(:,2), 'yymmdd'); 
    bl        = meancor_bl_dates.bl(idx); 


overwrite = 0; %overwrite current concatenated csv file? 

if overwrite == 1
    system(['rm ' csv_all]); 
    
    ixi=-1; 
for i=1:length(bl); 
    d1 = [d1_all(i,:)];
    d2 = [d2_all(i,:)];
    intdir = ['int_' d1 '_' d2]; 
    cd(intdir);
    
    
    % import csv
    csvi = importdata('meta_filt_topophase.flat.geo.csv'); 
    ni   = csvi.data; 
    nt   = csvi.textdata; 
    
    % print first line with column headers
    if i == 1
        fid = fopen(csv_all, 'at'); 
        fprintf(fid, 'id_no,idx,date1,date2,datenum1,datenum2,baseline\n'); 
    end
    
    % add data to csv file
    tifid = cell2mat(nt(2,1)); 
    ixi   = ixi+1; 
    fid = fopen(csv_all, 'at'); 
    fprintf(fid, [tifid ',']); 
    fprintf(fid, [num2str(ixi) ',']); 
    fprintf(fid, ['d' num2str(d1) ',' 'd' num2str(d2) ',' num2str(ni(1)) ',' ... 
        num2str(ni(2)) ',' num2str(ni(3)) '\n']); 
    fclose(fid); 
    
    % put geotiffs in folder
    system(['cp ' tifid '.tif ' geefol]); 
    
    cd ..
   
end
else 
    % do nothing 
end













