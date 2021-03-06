function gdalbuildvrt2iscexml(gdalvrtfile)

% Should be in WGS84 (use gdalwarp)
% Must be in folder containing gdal vrt file to translate
% gdalvrtfile e.g.: gdalvrtfile='outwgs84.dem.vrt'; 

% gdalvrtfile = 'outwgs84.dem.vrt'; 

[~, sztxt] = system(['more ' gdalvrtfile ' | grep raster']); 
    qidx   = strfind(sztxt, '"'); 
[~, dtxt]  = system(['more ' gdalvrtfile ' | grep GeoTransform']); 
    gtval  = regexp(dtxt,'[+-]?\d+\.?\d*e?[+-]?\d*', 'match')';

    
delta1          = cell2mat(gtval(2));
size1           = sztxt(qidx(1)+1:qidx(2)-1);
startingvalue1  = cell2mat(gtval(1));
endingvalue1    = num2str(str2double(startingvalue1)+(str2double(size1).*str2double(delta1)), '%e');
delta2          = cell2mat(gtval(6));
size2           = sztxt(qidx(3)+1:qidx(4)-1);
startingvalue2  = cell2mat(gtval(4));
endingvalue2    = num2str(str2double(startingvalue2)+(str2double(size2).*str2double(delta2)), '%e');
extra_file_vrt  = gdalvrtfile; 
file_name       = gdalvrtfile(1:end-4); 
width           = size1; 


filenamevrt = [file_name '.xml']; 
fid         = fopen(filenamevrt, 'wt');

fprintf(fid, '<imageFile> \n');
fprintf(fid, '    <property name="ISCE_VERSION"> \n'); 
fprintf(fid, '        <value>Release: 2.0.0_20170403, svn-2256, 20170403. Current: svn-Unknown.</value> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '    <property name="access_mode"> \n'); 
fprintf(fid, '        <value>read</value> \n'); 
fprintf(fid, '        <doc>Image access mode.</doc> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '    <property name="byte_order"> \n'); 
fprintf(fid, '        <value>l</value> \n'); 
fprintf(fid, '        <doc>Endianness of the image.</doc> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '    <component name="coordinate1"> \n'); 
fprintf(fid, '        <factorymodule>isceobj.Image</factorymodule> \n'); 
fprintf(fid, '        <factoryname>createCoordinate</factoryname> \n'); 
fprintf(fid, '        <doc>First coordinate of a 2D image (width).</doc> \n'); 
fprintf(fid, '        <property name="delta"> \n'); 
fprintf(fid, '            <value>%s</value> \n', delta1); 
fprintf(fid, '            <doc>Coordinate quantization.</doc> \n'); 
fprintf(fid, '        </property> \n'); 
fprintf(fid, '        <property name="endingvalue"> \n'); 
fprintf(fid, '            <value>%s</value> \n', endingvalue1); 
fprintf(fid, '            <doc>Starting value of the coordinate.</doc> \n'); 
fprintf(fid, '        </property> \n'); 
fprintf(fid, '        <property name="family"> \n'); 
fprintf(fid, '            <value>imagecoordinate</value> \n'); 
fprintf(fid, '            <doc>Instance family name</doc> \n'); 
fprintf(fid, '        </property> \n'); 
fprintf(fid, '        <property name="name"> \n'); 
fprintf(fid, '            <value>imagecoordinate_name</value> \n'); 
fprintf(fid, '            <doc>Instance name</doc> \n'); 
fprintf(fid, '        </property> \n'); 
fprintf(fid, '        <property name="size"> \n'); 
fprintf(fid, '            <value>%s</value> \n', size1); 
fprintf(fid, '            <doc>Coordinate size.</doc> \n'); 
fprintf(fid, '        </property> \n'); 
fprintf(fid, '        <property name="startingvalue"> \n'); 
fprintf(fid, '            <value>%s</value> \n', startingvalue1); 
fprintf(fid, '            <doc>Starting value of the coordinate.</doc> \n'); 
fprintf(fid, '        </property> \n'); 
fprintf(fid, '    </component> \n'); 
fprintf(fid, '    <component name="coordinate2"> \n'); 
fprintf(fid, '        <factorymodule>isceobj.Image</factorymodule> \n'); 
fprintf(fid, '        <factoryname>createCoordinate</factoryname> \n'); 
fprintf(fid, '        <doc>Second coordinate of a 2D image (length).</doc> \n'); 
fprintf(fid, '        <property name="delta"> \n'); 
fprintf(fid, '            <value>%s</value> \n', delta2); 
fprintf(fid, '            <doc>Coordinate quantization.</doc> \n'); 
fprintf(fid, '        </property> \n'); 
fprintf(fid, '        <property name="endingvalue"> \n'); 
fprintf(fid, '            <value>%s</value> \n', endingvalue2); 
fprintf(fid, '            <doc>Starting value of the coordinate.</doc> \n'); 
fprintf(fid, '        </property> \n'); 
fprintf(fid, '        <property name="family"> \n'); 
fprintf(fid, '            <value>imagecoordinate</value> \n'); 
fprintf(fid, '            <doc>Instance family name</doc> \n'); 
fprintf(fid, '        </property> \n'); 
fprintf(fid, '        <property name="name"> \n'); 
fprintf(fid, '            <value>imagecoordinate_name</value> \n'); 
fprintf(fid, '            <doc>Instance name</doc> \n'); 
fprintf(fid, '        </property> \n'); 
fprintf(fid, '        <property name="size"> \n'); 
fprintf(fid, '            <value>%s</value> \n', size2); 
fprintf(fid, '            <doc>Coordinate size.</doc> \n'); 
fprintf(fid, '        </property> \n'); 
fprintf(fid, '        <property name="startingvalue"> \n'); 
fprintf(fid, '            <value>%s</value> \n', startingvalue2); 
fprintf(fid, '            <doc>Starting value of the coordinate.</doc> \n'); 
fprintf(fid, '        </property> \n'); 
fprintf(fid, '    </component> \n'); 
fprintf(fid, '    <property name="data_type"> \n'); 
fprintf(fid, '        <value>short</value> \n'); 
fprintf(fid, '        <doc>Image data type.</doc> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '    <property name="extra_file_name"> \n'); 
fprintf(fid, '        <value>%s</value> \n', extra_file_vrt); 
fprintf(fid, '        <doc>For example name of vrt metadata.</doc> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '    <property name="family"> \n'); 
fprintf(fid, '        <value>demimage</value> \n'); 
fprintf(fid, '        <doc>Instance family name</doc> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '    <property name="file_name"> \n'); 
fprintf(fid, '        <value>%s</value> \n', file_name); 
fprintf(fid, '        <doc>Name of the image file.</doc> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '    <property name="image_type"> \n'); 
fprintf(fid, '        <value>dem</value> \n'); 
fprintf(fid, '        <doc>Image type used for displaying.</doc> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '    <property name="length"> \n'); 
fprintf(fid, '        <value>%s</value> \n', size2); 
fprintf(fid, '        <doc>Image length</doc> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '    <property name="name"> \n'); 
fprintf(fid, '        <value>demimage_name</value> \n'); 
fprintf(fid, '        <doc>Instance name</doc> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '    <property name="number_bands"> \n'); 
fprintf(fid, '        <value>1</value> \n'); 
fprintf(fid, '        <doc>Number of image bands.</doc> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '    <property name="reference"> \n'); 
fprintf(fid, '        <value>WGS84</value> \n'); 
fprintf(fid, '        <doc>Geodetic datum</doc> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '    <property name="scheme"> \n'); 
fprintf(fid, '        <value>BIP</value> \n'); 
fprintf(fid, '        <doc>Interleaving scheme of the image.</doc> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '    <property name="width"> \n'); 
fprintf(fid, '        <value>%s</value> \n', width); 
fprintf(fid, '        <doc>Image width</doc> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '   <property name="xmax"> \n'); 
fprintf(fid, '        <value>%s</value> \n', endingvalue1); 
fprintf(fid, '        <doc>Maximum range value</doc> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '    <property name="xmin"> \n'); 
fprintf(fid, '        <value>%s</value> \n', startingvalue1); 
fprintf(fid, '        <doc>Minimum range value</doc> \n'); 
fprintf(fid, '    </property> \n'); 
fprintf(fid, '</imageFile> \n'); 

fclose(fid);


