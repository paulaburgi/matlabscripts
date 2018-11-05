
dates = [20070713;20070828;20071013;20071128;20080113;20080228;20080414;...
         20080530;20080715;20080830;20081130;20090115;20090302;20090718;...
         20091018;20100118;20100420;20100605;20100721;20101206;20110308];
nd=length(dates);
d = num2str(dates); 

geefold = 'for_gee/'; 
% download data
for i=1:nd
    if ~exist([geefold d(i,:) '.tif'], 'file')
    string=['wget https://mesonet.agron.iastate.edu/request/gis/n0r2gtiff.php' ...
            '?dstr=' d(i,:) '0635 --output-document=' d(i,:) '_nexrad.zip'];
    system(string);
    system(['unzip ' d(i,:) '_nexrad.zip']);
    system(['gdal_translate -projwin -124 45 -122 43 n0r_' d(i,:) '0635.tif ' geefold d(i,:) '.tif']);
    end
end


% make CSV metadata
fid = fopen([geefold 'meta_nexrad_v1.csv'], 'wt'); 
fprintf(fid, 'id_no,idx,dt\n'); 
fclose(fid); 
ix = 0; 
for i=1:nd
    fid = fopen([geefold 'meta_nexrad_v1.csv'], 'a'); 
    fprintf(fid, [d(i,:) ',']); 
    fprintf(fid, [num2str(ix) ',']); 
    fprintf(fid, ['d' d(i,:) '\n']); 
    fclose(fid); 
    ix = ix+1; 
end


% for i=1:nd
%     string=['wget https://mesonet.agron.iastate.edu/archive/data/' ... 
%             d(i,1:4) '/' d(i,5:6) '/' d(i,7:8) '/GIS/uscomp/n0r_' d(i,:) '0650.tif' ...
%             ' --output-document=' d(i,:) '_nexrad.zip'];
%     system(string);
%     %system(['unzip ' d(i,:) '_nexrad.zip']);
%     system(['gdal_translate -projwin -124 45 -122 43 n0r_' ...
%             d(i,:) '0650.tif  ' d(i,:) '.tif']);
% end

%https://mesonet.agron.iastate.edu/archive/data/%Y/%m/%d/GIS/uscomp/n0r_%Y%m%d%H%M.png

% string=['wget https://mesonet.agron.iastate.edu/archive/data/2007/07/13/GIS/uscomp/n0r_200707130650.png ' ...
%         '--output-document=test_nexrad.zip'];