t  = importdata('/data/pmb229/isce/p222f870/data/analysis/geotiff_gee/meta_all.csv'); 
td = cell2mat(t.textdata(2:end,3:4));  
nints = length(td); 
intdirs = [repmat('int_', nints, 1) td(:,2:7) repmat('_', nints, 1) td(:,9:end)];
vrtdir = 'aligned_insar'; 
%vrtdir = 'VRT'; 

% put all unw and cor vrts in 1 folder
%if length(dir([vrtdir '/*vrt'])) ~= length(intdirs).*2
    for i=1:length(intdirs)
        intdir = intdirs(i,:); 
        dates = ['20' intdir(5:10) '_' '20' intdir(12:end)];
        

        system(['cp ' intdir '/filt_topophase.unw.geo.vrt ' vrtdir '/' ...
                dates '_unw.vrt']); 
        system(['cp ' intdir '/phsig.cor.geo.vrt ' vrtdir '/' ...
                dates '_cor.vrt']); 

        system(['sed -i "s|filt_topophase.unw.geo|/data/pmb229/isce/p222f870/mostcombos/' ...
            intdir '/filt_topophase.unw.geo.vrt|g" ' vrtdir '/' dates '_unw.vrt']);
        system(['sed -i "s|phsig.cor.geo|/data/pmb229/isce/p222f870/mostcombos/' ...
            intdir '/phsig.cor.geo.vrt|g" ' vrtdir '/' dates '_cor.vrt']);
        
       % system(['ln -sf ../' intdir '/filt_topophase.unw.geo ' vrtdir '/'dates '_]); 
    end
%end


% put all isce logs in 1 folder
% logfol = 'iscelogs'; 
% for i=1:length(intdirs)
%         intdir = intdirs(i,:); 
%         dates = ['20' intdir(5:10) '_' '20' intdir(12:end)]; 
% 
%         system(['cp ' intdir '/isce.log ' logfol '/' ...
%                 dates '.isce.log']); 
% end
    


















