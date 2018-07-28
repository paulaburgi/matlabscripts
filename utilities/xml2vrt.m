% xml2vrt for filt_topophase.flat.geo.int

% oregon
pf_fol  = '/data/pmb229/isce/p222f870/ALOS2_data/';
cd(pf_fol); 
intdirs = dir; 
intdirs = {intdirs.name}; 

for i = 1:length(intdirs)
    intdir = cell2mat(intdirs(i)); 
    if contains(intdir, 'int_')
        cd(intdir)

% for some reason, matlab isn't reading in xml file, must convert to txtfile
xmlfile = 'filt_fine_masked.6alks_3rlks.int.geo.xml'; 
% xmlfile = 'filt_topophase.flat.geo.xml';
txtfile = [xmlfile(1:end-3) 'txt']; 
t   = system(['cp ' xmlfile ' ' txtfile]);
xml = importdata(txtfile); 
t2  = system(['rm ' txtfile]);

sz = strfind(xml, '"size"'); 

% sizes
[~,bltxt] = system(['more ' xmlfile ' | grep size -A 1']);
digidx = regexp(bltxt, '\d');
consec = diff(digidx)==1; 
endidx = find(consec == 0); 
size1  = bltxt(digidx(1:endidx));
size2  = bltxt(digidx(endidx+1:end));

% starting values
[~,bltxt] = system(['more ' xmlfile ' | grep startingvalue -A 1']);
digidx = regexp(bltxt, '[-\d.]');
numstr = bltxt(digidx); 
ddfind = strfind(numstr, '--'); 
digidx = digidx([1:ddfind-1 ddfind+2:end]);
consec = diff(digidx)==1;
endidx = find(consec == 0); 
start1  = bltxt(digidx(1:endidx(1)));
start2  = bltxt(digidx(endidx+1:end));

% delta
[~,bltxt] = system(['more ' xmlfile ' | grep delta -A 1']);
digidx = regexp(bltxt, '[-\d.]');
numstr = bltxt(digidx); 
ddfind = strfind(numstr, '--'); 
digidx = digidx([1:ddfind-1 ddfind+2:end]);
consec = diff(digidx)==1;
endidx = find(consec == 0); 
delta1  = bltxt(digidx(1:endidx(1)));
delta2  = bltxt(digidx(endidx+1:end));

% file name
[~,bltxt] = system(['more ' xmlfile ' | grep file_name -A 1']);
valfind1  = strfind(bltxt, '<value>');
valfind2  = strfind(bltxt, '</value>');
lv        = length(valfind1);
if lv == 1; 
    idx = 1; 
else
    idx = 2; 
end
file_name = bltxt(valfind1(idx)+7:valfind2(idx)-1);

% line offset
pixeloff = '8'; 
lineoff  = num2str(str2num(pixeloff).*str2num(size1));


% write file 
filenamevrt = [file_name '.vrt']; 
fid         = fopen(filenamevrt, 'wt');

fprintf(fid, ['<VRTDataset rasterXSize="' size1 '" rasterYSize="' size2 '">\n']); 
fprintf(fid, ['    <SRS>EPSG:4326</SRS>\n']); 
fprintf(fid, ['    <GeoTransform>' start1 ', ' delta1 ', 0.0, ' start2 ', 0.0, ' delta2 '</GeoTransform>\n']); 
fprintf(fid, ['    <VRTRasterBand band="1" dataType="CFloat32" subClass="VRTRawRasterBand">\n']); 
fprintf(fid, ['        <SourceFilename relativeToVRT="1">' file_name '</SourceFilename>\n']); 
fprintf(fid, ['        <ByteOrder>LSB</ByteOrder>\n']); 
fprintf(fid, ['        <ImageOffset>0</ImageOffset>\n']); 
fprintf(fid, ['        <PixelOffset>' pixeloff '</PixelOffset>\n']); 
fprintf(fid, ['        <LineOffset>' lineoff '</LineOffset>\n']); 
fprintf(fid, ['    </VRTRasterBand>\n']); 
fprintf(fid, ['</VRTDataset>']); 

fclose(fid);
cd ..

    else 
        %pass
    end
    
end





