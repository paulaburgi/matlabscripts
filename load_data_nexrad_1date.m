
% dates = [20070713;20070828;20071013;20071128;20080113;20080228;20080414;...
%          20080530;20080715;20080830;20081130;20090115;20090302;20090718;...
%          20091018;20100118;20100420;20100605;20100721;20101206;20110308];
% nd=length(dates);
% d = num2str(dates); 

cd /data/pmb229/other/nexrad/nexrad_080414

mins = [0:5:55]'; 
hrs  = [0:6];
hm = {};
for j=1:length(hrs); 
    hri = num2str(hrs(j)); 
    if length(hri) == 1
        hri = ['0' hri]; 
    end
    for k = 1:length(mins)
        mini = num2str(mins(k)); 
        if length(mini) == 1
            mini = ['0' mini]; 
        end
        idx = (j-1)*length(mins) + k;
        hm{idx,1} = [hri mini]; 
    end
end
maxtime = '0635'; 
idx     = find(strcmp(hm, maxtime) == 1); 
hm      = cell2mat({hm{1:idx,1}}'); 
ntimes = length(hm); 

% date
d = '20080414'; 

geefold = 'for_gee/'; 
% download data
for i=1:ntimes
    itime = hm(i,:); 
    if ~exist([geefold d '_' itime '.tif'], 'file')
        string=['wget https://mesonet.agron.iastate.edu/request/gis/n0r2gtiff.php' ...
                '?dstr=' d itime ' --output-document=' d '_' itime '_nexrad.zip'];
        system(string);
        system(['unzip ' d '_' itime '_nexrad.zip']);
        system(['gdal_translate -projwin -124 45 -122 43 n0r_' d itime '.tif ' geefold d '_' itime '.tif']);
    end
end


% make CSV metadata
fid = fopen([geefold 'meta_nexrad_' d '_v1.csv'], 'wt'); 
fprintf(fid, 'id_no,idx,dt\n'); 
fclose(fid); 
ix = 0; 
for i=1:ntimes
    fid = fopen([geefold 'meta_nexrad_' d '_v1.csv'], 'a'); 
    fprintf(fid, [d '_' hm(i,:) ',']); 
    fprintf(fid, [num2str(ix) ',']); 
    fprintf(fid, ['d' d hm(i,:) '\n']); 
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