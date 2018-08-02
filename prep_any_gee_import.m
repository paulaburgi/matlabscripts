% prep_any_gee_import

clear 
close all

f = 'all_stitched.dem.geo'; 

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
    elseif contains(f, 'all_stitched')
        midx      = find(fg == -9988 | fg == 0); 
    end
    fg(midx)  = ndval;
    
    % write tif file 
    fh      = strrep(f, '.', '_'); 
    tifname = ['gee_' fh '.tif']; 
    geotiffwrite(tifname, flipud(fg), info.RefMatrix);
    
    











