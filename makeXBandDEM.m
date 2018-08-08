dtdir = dir('*dt2'); 
dtdir = {dtdir.name}; 


for i=1:length(dtdir); 
    dti = cell2mat(dtdir(i)); 
    demi = [dti '.dem'];
    if ~exist(demi, 'file'); 
        system(['gdal_translate -of ENVI ' dti ' ' demi]); 
        system(['gdalbuildvrt ' demi '.vrt ' demi]); 
    end
end


system(['gdalbuildvrt stitchedDEM.vrt *vrt']); 
system(['gdal_translate -of ENVI stitchedDEM.vrt stitched.dem']); 

