% OregonLidar2isceDEM.m

cdir  = pwd; 
% dldir = [pwd '/downlooked2/'];
% hdir  = 'Bare_Earth/';
% prfx  = 'be'; 
dldir = [pwd '/hhdownlooked/']; 
hdir  = 'Highest_Hit/';
prfx  = 'hh'; 

extfol = [cdir '/extrafiles']; %sta.adf is a single image's stats. 

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
    indir = [prfx Ldiri(5:9) lower(Ldiri(10)) Ldiri(11)]; 
    indir2= [indir(1:2) '_' indir(3:end)]; 
     
    cd(Ldiri)
%     ddir = dir('2015*');
    ddir = dir('20*');
    ddir = {ddir.name};
    zipc = contains(ddir, 'zip'); 
        for m=1:length(ddir)
            if zipc(m) == 1
                zipcm = cell2mat(ddir(m)); 
                if contains(zipcm, ' '); 
                    sidx = strfind(zipcm, ' '); 
                    zipcm = [zipcm(1:sidx-1) '\' zipcm(sidx:end)];
                end
                system(['unzip ' zipcm]);
                system(['rm ' zipcm]); 
            end
        end
    ddir = dir('20*');
    ddir = {ddir.name};

    for j=1:length(ddir)
        sdt   = cell2mat(ddir(j));
        cd(sdt); 
        sdtdir = dir('*zip'); 
        sdtdir = {sdtdir.name};
        if ~isempty(sdtdir)
            for n=1:length(sdtdir)
                sdtdirn = cell2mat(sdtdir(n)); 
                system(['unzip ' sdtdirn]); 
                system(['rm ' sdtdirn]); 
            end
        end
        if exist(hdir, 'file')
            cd(hdir);
        else
            cd(strrep(hdir, '_', ''));
        end
        tfile = '*.img'; 
        imgs  = length(dir(tfile));
        if exist(indir, 'file') || imgs > 0 || exist(indir2, 'file')
            if exist(indir, 'file')
                cd(indir);
                tfile = 'w001001.adf';
            elseif exist(indir2, 'file')
                cd(indir2); 
                tfile = 'w001001.adf'; 
            end
            if ~exist('prj.adf', 'file'); 
                system(['cp ' extfol ' .']); 
                system(['cp ' extfol '/* .']); 
                % ADD MORE FILES HERE SO YOU STOP GETTING ERROR
            end
            if ~exist('outwgs84.dem', 'file'); 
                %system('rm out*'); 
                system(['gdal_translate -ot Int16 -of ENVI ' tfile ' out.dem']); 
                system(['gdalbuildvrt out.dem.vrt out.dem']); 
                system(['gdalwarp -t_srs EPSG:4326 -of ENVI -r bilinear out.dem.vrt outwgs84.dem']); % if it fails, cp out.dem.vrt from be...
                system(['gdalbuildvrt outwgs84.dem.vrt outwgs84.dem']); 
                gdalbuildvrt2iscexml('outwgs84.dem.vrt');
            end
            if ~exist('outwgs84m.dem', 'file');
                %system('rm outwgs84m.d*'); 
                system(['imageMath.py -e=''a*0.3048'' -o outwgs84m1.dem  --a=outwgs84.dem']); % convert ft to m
                system(['gdal_translate -ot Int16 -of ENVI outwgs84m1.dem outwgs84m.dem']); 
                system(['gdalbuildvrt outwgs84m.dem.vrt outwgs84m.dem']); 
                gdalbuildvrt2iscexml('outwgs84m.dem.vrt');
            end
            foldnm  = strrep(sdt, ' ', '_'); 
            demfile = [indir '_' foldnm '.dem'];
            if ~exist([dldir demfile], 'file'); 
                system(['looks.py -i outwgs84m.dem -o ' demfile ' -r 10 -a 10']);
                system(['mv ' prfx '* ' dldir]); 
            end
        end 
            
        if exist([indir 'a/'], 'file') || exist([indir 'b/'], 'file') || exist([indir '1/'], 'file') || exist([indir '2/'], 'file'); 
            clear abdir
            abdir  = dir([indir '*']); 
                z  = cell2mat({abdir.isdir}); 
            abdir  = abdir(z); 
            abdir  = {abdir.name};
            abdirl = length(abdir);
            tfile = 'w001001.adf';
            for k = 1:abdirl
                indira = cell2mat(abdir(k)); 
                cd(indira); 
                if ~exist('outwgs84.dem', 'file'); 
                    %system('rm out*'); 
                    system(['gdal_translate -ot Int16 -of ENVI ' tfile ' out.dem']); 
                    system(['gdalbuildvrt out.dem.vrt out.dem']); 
                    system(['gdalwarp -t_srs EPSG:4326 -of ENVI -r bilinear out.dem.vrt outwgs84.dem']); 
                    system(['gdalbuildvrt outwgs84.dem.vrt outwgs84.dem']); 
                    gdalbuildvrt2iscexml('outwgs84.dem.vrt');
                end
                if ~exist('outwgs84m.dem', 'file');
                    %system('rm outwgs84m.d*'); 
                    system(['imageMath.py -e=''a*0.3048'' -o outwgs84m1.dem  --a=outwgs84.dem']); % convert ft to m
                    system(['gdal_translate -ot Int16 -of ENVI outwgs84m1.dem outwgs84m.dem']); 
                    system(['gdalbuildvrt outwgs84m.dem.vrt outwgs84m.dem']); 
                    gdalbuildvrt2iscexml('outwgs84m.dem.vrt');
                end
                foldnm  = strrep(sdt, ' ', '_'); 
                demfile = [indir '_' foldnm '_' num2str(k) '.dem'];
                if ~exist([dldir demfile], 'file'); 
                    system(['looks.py -i outwgs84m.dem -o ' demfile ' -r 10 -a 10']);
                    system(['mv ' prfx '* ' dldir]); 
                end
                cd ..
            end
        end
            
            
            
            
        cd([cdir '/' Ldiri]);
    end  
    cd(cdir);
end






% % stitch the dem together 
    cd(dldir);
    system('rm all*'); 
    system('gdalbuildvrt -srcnodata -9988 all.dem.vrt *vrt'); 
    system('gdal_translate -of ENVI -a_nodata -9988 all.dem.vrt all_stitched.dem'); 
    system('gdalbuildvrt -srcnodata -9988 all_stitched.dem.vrt all_stitched.dem');
    gdalbuildvrt2iscexml('all_stitched.dem.vrt');

    
    
    
    
% % mask the no data values as nan 
    
    filename = 'all_stitched.dem';
    vrtfilename = [filename '.vrt']; 
    [~, sztxt]  = system(['more ' vrtfilename ' | grep raster']); 
    qidx        = strfind(sztxt, '"'); 
    nx          = str2double(sztxt(qidx(1)+1:qidx(2)-1));
    ny          = str2double(sztxt(qidx(3)+1:qidx(4)-1));
    fid         = fopen(filename,'r','native');
    [rmg,count] = fread(fid,[nx,ny],'int16');
    rmg2 = rmg;
    idx = find(rmg < 0);
    rmg2(idx) = 1e6;
    fid = fopen('all_stitched_masked.dem', 'wb');
    fwrite(fid, rmg2, 'int16');
    fclose(fid);
    system(['cp all_stitched.dem.vrt ' filename '_masked.dem.vrt']); 
    system(['cp all_stitched.dem.xml ' filename '_masked.dem.xml']); 
    system('sed -i s/all_stitched/all_stitched_masked/g all_stitched_masked.dem.xml'); 
        % above replaces any occurence of "all_stitched" with"all_stitched_masked"
    
    
    
    
    
% gdal_translate -ot Int16 -of ENVI w001001.adf out.dem
% gdalbuildvrt out.dem.vrt out.dem
% gdalwarp -t_srs EPSG:4326 -of ENVI -r bilinear out.dem.vrt outwgs84.dem
% gdalbuildvrt outwgs84.dem.vrt outwgs84.dem
%     
%     gdalbuildvrt2iscexml('outwgs84.dem.vrt');
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    