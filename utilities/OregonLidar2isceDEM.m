% OregonLidar2isceDEM.m

cdir = pwd; 
dldir = [pwd '/downlooked2/'];


% unzip any zipped lidar folders
zdir = dir('LDQ*zip'); 
zdir = {zdir.name};

if ~isempty(zdir) 
    for i = 1:length(zdir)
        zdiri = cell2mat(zdir(i)); 
        system(['unzip ' zdiri]); 
        system(['rm ' zdiri]);
    end
end

Ldir = dir('LDQ*');
Ldir = {Ldir.name};

Ll = length(Ldir); 

for i = 1:Ll 
    Ldiri = cell2mat(Ldir(i)); 
    indir = ['be' Ldiri(5:9) lower(Ldiri(10)) Ldiri(11)]; 
     
    cd(Ldiri)
%     ddir = dir('2015*');
    ddir = dir('20*');
    ddir = {ddir.name};

    for j=1:length(ddir)
        sdt   = cell2mat(ddir(j));
        cd([sdt '/Bare_Earth/']); 
        tfile = '*.img'; 
        imgs  = length(dir(tfile));
        if exist(indir, 'file') || imgs > 0 ; 
            if exist(indir, 'file')
                cd(indir);
                tfile = 'w001001.adf';
            end
            if ~exist('outwgs84.dem', 'file'); 
                system(['gdal_translate -ot Int16 -of ENVI -a_nodata -32768 ' tfile ' out.dem']); 
                system(['gdalbuildvrt -srcnodata -32768 out.dem.vrt out.dem']); 
                system(['gdalwarp -t_srs EPSG:4326 -of ENVI -r bilinear out.dem.vrt outwgs84.dem']); 
                system(['gdalbuildvrt outwgs84.dem.vrt outwgs84.dem']); 
                gdalbuildvrt2iscexml('outwgs84.dem.vrt');
            end
            foldnm  = strrep(sdt, ' ', '_'); 
            demfile = [indir '_' foldnm '.dem'];
            if ~exist([dldir demfile], 'file'); 
                system(['looks.py -i outwgs84.dem -o ' demfile ' -r 10 -a 10']);
                system(['mv be* ' dldir]); 
            end
        end 
            
        if exist([indir 'a/'], 'file') || exist([indir 'b/'], 'file') ; 
            clear abdir
            abdir  = dir([indir '*']); 
                z  = cell2mat({abdir.isdir}); 
            abdir  = abdir(z); 
            abdir  = {abdir.name};
            abdirl = length(abdir);
            for k = 1:abdirl
                indira = cell2mat(abdir(k)); 
                cd(indira); 
                if ~exist('outwgs84.dem', 'file'); 
                    system(['gdal_translate -ot Int16 -of ENVI -a_nodata -32768 ' tfile ' out.dem']); 
                    system(['gdalbuildvrt -srcnodata -32768 out.dem.vrt out.dem']); 
                    system(['gdalwarp -t_srs EPSG:4326 -of ENVI -r bilinear out.dem.vrt outwgs84.dem']); 
                    system(['gdalbuildvrt outwgs84.dem.vrt outwgs84.dem']); 
                    gdalbuildvrt2iscexml('outwgs84.dem.vrt');
                end
                foldnm  = strrep(sdt, ' ', '_'); 
                demfile = [indir '_' foldnm '_' num2str(k) '.dem'];
                if ~exist([dldir demfile], 'file'); 
                    system(['looks.py -i outwgs84.dem -o ' demfile ' -r 10 -a 10']);
                    system(['mv be* ' dldir]); 
                end
                cd ..
            end
        end    
            
            
            
            
        cd([cdir '/' Ldiri]);
    end  
    cd(cdir);
end




% % stitch the dem together
%     cd(dldir); 
%     system('rm all*'); 
%     system('gdalbuildvrt -srcnodata -32768 all.dem.vrt *vrt'); 
%     system('gdal_translate -of ENVI -a_nodata -32768 all.dem.vrt all_stitched.dem'); 
%     system('gdalbuildvrt -srcnodata -32768 all_stitched.dem.vrt all_stitched.dem');
%     gdalbuildvrt2iscexml('all_stitched.dem.vrt');

    
    
    
    
% % mask the no data values as nan 
    
%     filename = 'all_stitched.dem';
%     vrtfilename = [filename '.vrt']; 
%     [~, sztxt]  = system(['more ' vrtfilename ' | grep raster']); 
%     qidx        = strfind(sztxt, '"'); 
%     nx          = str2double(sztxt(qidx(1)+1:qidx(2)-1));
%     ny          = str2double(sztxt(qidx(3)+1:qidx(4)-1));
%     fid         = fopen(filename,'r','native');
%     [rmg,count] = fread(fid,[nx,ny],'int16');
%     rmg2 = rmg;
%     idx = find(rmg < 0);
%     rmg2(idx) = NaN;
%     fid = fopen('all_stitched_masked.dem', 'wb');
%     fwrite(fid, rmg2, 'int16');
%     fclose(fid);
%     system('cp all_stitched.dem.vrt all_stitched_masked.dem.vrt'); 
%     system('cp all_stitched.dem.xml all_stitched_masked.dem.xml'); 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    