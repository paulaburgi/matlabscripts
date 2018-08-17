% prep_any_gee_import

clear 
close all

%f = 'demdiff3.geo'; 
f = 'all_stitched.dem'; 
%f = 'all_stitched_masked.dem'; 
%f = 'demLat_N43_N45_Lon_W125_W121.dem.wgs84'; 
%f = 'z.rdr.geo'; 

% write a tif file, to get ref frame 
    if ~exist([f '.tif'], 'file'); 
        system(['gdal_translate ' f '.vrt ' f '.tif']); 
    end
    info = geotiffinfo([f '.tif']); 
    fg = geotiffread([f '.tif']); 
    

% identify no data areas 
    ndval     = -9999;
    if contains(f, 'demdiff'); 
        midx1     = find(fg == 0); 
        fg(midx1) = ndval; 
        fg = fg+24; 
        midx      = find(fg > 200 | fg < -200); 
    elseif contains(f, 'all_stitched_masked')
        midx      = find(fg > 3e4);
    elseif contains(f, 'all_stitched')
        midx      = find(fg == -9988); 
    elseif contains(f, 'demLat_') || contains(f, 'z.rdr')
        midx      = find(fg < 0); 
    end
    fg(midx)  = ndval;
    
    % write tif file 
    fh      = strrep(f, '.', '_'); 
    tifname = ['for_gee/gee_' fh '.tif']; 
    geotiffwrite(tifname, fg, info.RefMatrix);
    
    











