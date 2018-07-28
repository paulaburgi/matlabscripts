%gdalbuildvrt2iscevrt.m
%actually I don't think you need this, the gdal vrt works. 

iscevrtfile = 'be44123a4.dem.vrt';
gdalvrtfile = 'outwgs84.dem.vrt'; 

[~, sztxt] = system(['more ' gdalvrtfile ' | grep raster']); 
    qidx   = strfind(sztxt, '"'); 
[~, dtxt]  = system(['more ' gdalvrtfile ' | grep GeoTransform']); 
    gtval  = regexp(dtxt,'[+-]?\d+\.?\d*e?[+-]?\d*', 'match')';



size1       = sztxt(qidx(1)+1:qidx(2)-1);
size2       = sztxt(qidx(3)+1:qidx(4)-1);
startval1   = cell2mat(gtval(1));
startval2   = cell2mat(gtval(4));
delta1      = cell2mat(gtval(2));
delta2      = cell2mat(gtval(6));
demname     = iscevrtfile(1:end-4); 
pixoffset   = num2str(4); 
lineoffset  = num2str(str2double(pixoffset).*str2double(size1)); %pixeloffset*rasterXSize


fid         = fopen(iscevrtfile, 'wt');

fprintf(fid, '<VRTDataset rasterXSize="%s" rasterYSize="%s"> \n', size1, size2);
fprintf(fid, '    <SRS>EPSG:4326</SRS> \n');
fprintf(fid, '    <GeoTransform>%s, %s, 0.0, %s, 0.0, %s</GeoTransform> \n', startval1, delta1, startval2, delta2);
fprintf(fid, '    <VRTRasterBand band="1" dataType="Int16" subClass="VRTRawRasterBand"> \n');
fprintf(fid, '        <SourceFilename relativeToVRT="1">%s</SourceFilename> \n', demname);
fprintf(fid, '        <ByteOrder>LSB</ByteOrder> \n');
fprintf(fid, '        <ImageOffset>0</ImageOffset> \n');
fprintf(fid, '        <PixelOffset>%s</PixelOffset> \n', pixoffset);
fprintf(fid, '        <LineOffset>%s</LineOffset> \n', lineoffset);
fprintf(fid, '    </VRTRasterBand> \n');
fprintf(fid, '</VRTDataset> \n');

fclose(fid);
