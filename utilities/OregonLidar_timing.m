% OregonLidar2isceDEM.m

cd('/data/pmb229/other/OregonLidarDEM');
cdir = pwd; 
cd('downlooked2/'); 
tdir = '../DEMtiming/'; 


% unzip any zipped lidar folders
bdir = dir('b*dem'); 
bdir = {bdir.name};

for i=1:length(bdir)
    clear rmg
    bfile = cell2mat(bdir(i)); 
    yfile       = ['../DEMtiming/t_' bfile]; 
    ydlfile     = ['../DEMtiming/tdl_' bfile]; 
    
    
    if ~exist(ydlfile, 'file')
        % load DEM
        hidx  = strfind(bfile, '_'); 
        yr    = str2double(bfile(hidx(1)+1:hidx(2)-1));  
        vrtfilename = [bfile '.vrt']; 
        [~, sztxt]  = system(['more ' vrtfilename ' | grep raster']); 
        qidx        = strfind(sztxt, '"'); 
        nx          = str2double(sztxt(qidx(1)+1:qidx(2)-1));
        ny          = str2double(sztxt(qidx(3)+1:qidx(4)-1));
        fid         = fopen(bfile,'r','native');
        [rmg,count] = fread(fid,[nx,ny],'int16');

        % change values to year acquired and write file
        idx         = find(rmg ~= -9988); 
        rmg(idx)    = yr;
        fid = fopen(yfile, 'wb');
        fwrite(fid, rmg, 'int16');
        fclose(fid);

        % copy and change metadata
        system(['cp ' bfile '.vrt ' yfile '.vrt']);
        system(['cp ' bfile '.xml ' yfile '.xml']);
        system(['sed -i s/' bfile '/' yfile(14:end) '/g ' yfile '.xml']); 
        system(['sed -i s/' bfile '/' yfile(14:end) '/g ' yfile '.vrt']); 
        system(['chmod 777 ' yfile '*']);
        system(['looks.py -i ' yfile ' -o ' ydlfile  ' -r 10 -a 10']);
        
        % remove non-downlooked files
        system(['rm ' yfile '*']); 
    end
    
end
  


% % stitch the dem together
    cd(tdir)
    if ~exist('allyear_stitched.dem.xml'); 
        system('gdalbuildvrt -srcnodata -9988 allyear.dem.vrt tdl*vrt'); 
        system('gdal_translate -of ENVI -a_nodata -9988 allyear.dem.vrt allyear_stitched.dem'); 
        system('gdalbuildvrt -srcnodata -9988 allyear_stitched.dem.vrt allyear_stitched.dem');
        gdalbuildvrt2iscexml('allyear_stitched.dem.vrt');
    end
    
% % mask the no data values as nan 
    
    if ~exist('allyear_stitched_masked.dem.xml'); 
        filename = 'allyear_stitched.dem';
        vrtfilename = [filename '.vrt']; 
        [~, sztxt]  = system(['more ' vrtfilename ' | grep raster']); 
        qidx        = strfind(sztxt, '"'); 
        nx          = str2double(sztxt(qidx(1)+1:qidx(2)-1));
        ny          = str2double(sztxt(qidx(3)+1:qidx(4)-1));
        fid         = fopen(filename,'r','native');
        [rmg,count] = fread(fid,[nx,ny],'int16');
        rmg2 = rmg;
        idx = find(rmg < 2000);
        rmg2(idx) = -9988;
        fid = fopen('allyear_stitched_masked.dem', 'wb');
        fwrite(fid, rmg2, 'int16');
        fclose(fid);
        system('cp allyear_stitched.dem.vrt allyear_stitched_masked.dem.vrt'); 
        system('cp allyear_stitched.dem.xml allyear_stitched_masked.dem.xml'); 
        system(['sed -i s/allyear_stitched/allyear_stitched_masked/g allyear_stitched.dem.xml']); 
        system(['sed -i s/allyear_stitched/allyear_stitched_masked/g allyear_stitched.dem.vrt']);
    end

% put in gee
    f = 'allyear_stitched_masked.dem'; 
    if ~exist('for_gee', 'file'); 
        system('mkdir for_gee'); 
    end
    if ~exist('for_gee/*.tif', 'file'); 
        system(['gdal_translate ' f '.vrt ' f '.tif']); 
        info = geotiffinfo([f '.tif']); 
    fg = geotiffread([f '.tif']); 
    ndval     = -9999;
    midx      = find(fg == -9988); 
    fg(midx)  = ndval;
    fh      = strrep(f, '.', '_'); 
    tifname = ['for_gee/gee_' fh '.tif']; 
    geotiffwrite(tifname, fg, info.RefMatrix);
%     system(['geebam upload -u paula.burgi@gmail.com --source ' ... 
%         cdir '/DEMtiming/for_gee ' ...
%         '--dest users/pmb229/cascadia/OregonLidarTiming --nodata -9999']); 
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    