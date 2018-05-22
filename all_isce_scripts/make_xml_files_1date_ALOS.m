% make .xml files for ALOS data, to get baseline data
% ******make sure you're in folder data/baselines******

clear

% data folder
  % oregon 
    % datafol = '/data/pmb229/isce/p222f870/data/'; 
    % baselinefol = [datafol 'baselines/']; 
  % sumatra
    datafol = '/data/pmb229/isce/p446f7190_sumatra/data/'; 
    baselinefol = [datafol 'baselines/']; 

    cd(baselinefol)
    s = dir('../'); 


% parameters
    filter_strength = 0.3;  % 0.0 - 1.0
    unwrap = 'False'; % 'True' or 'False'


% find file/folders with digits in file/folder name
    f=[]; 
    for i=1:length(s)
        TF = isstrprop(s(i).name,'digit'); 
        f = [f; sum(TF)];
    end

% find indices of files slc folders 
    slcidx = find(f==6);
    nslc   = length(slcidx);

% concatenate all slc folder names
    fids=[]; 
    for i=1:nslc
        sn = s(slcidx(i)).name; 
        fids=[fids; sn];
    end

% first date, baselines will be relative to
    s1 = fids(1,:); 

% make xml files
    for i = 2:nslc
        s2 = fids(i,:); 
        intfol = sprintf('int_%s_%s', s1, s2); 
        datafol_s1 = sprintf('%s/%s', datafol, s1); 
        datafol_s2 = sprintf('%s/%s', datafol, s2); 
        s11 = dir(datafol_s1); 
        s21 = dir(datafol_s2); 
        s12 = {s11.name}; 
        s22 = {s21.name}; 
        s13 = cell2mat(s12); 
        s23 = cell2mat(s22); 
        IMG1 = cell2mat(s12(contains(s12, 'IMG-HH'))); 
        IMG2 = cell2mat(s22(contains(s22, 'IMG-HH'))); 
        LED1 = cell2mat(s12(contains(s12, 'LED'))); 
        LED2 = cell2mat(s22(contains(s22, 'LED'))); 

        if ~exist(intfol, 'file'); 
            system(['mkdir ' intfol]);
        end
        cd(intfol); 

        %link data files
        system(sprintf('ln -sf %s/%s/%s %s', datafol, s1, IMG1, IMG1)); 
        system(sprintf('ln -sf %s/%s/%s %s', datafol, s2, IMG2, IMG2)); 
        system(sprintf('ln -sf %s/%s/%s %s', datafol, s1, LED1, LED1)); 
        system(sprintf('ln -sf %s/%s/%s %s', datafol, s2, LED2, LED2));

        % write xml
        fid=fopen(sprintf('%s.xml', intfol), 'wt'); 
        fprintf(fid, '<insarApp>\n'); 
        fprintf(fid, '	<component name="insar">\n\n'); 
        fprintf(fid, '       <property name="Sensor Name">\n'); 
        fprintf(fid, '           <value>ALOS</value>\n'); 
        fprintf(fid, '       </property>\n\n'); 
        fprintf(fid, '<!--       uncomment this "doppler method" if you get an error related to it-->\n'); 
        fprintf(fid, '<!--       <property name="doppler method">useDEFAULT</property>-->\n\n'); 
        fprintf(fid, '<!--   you can use either posting or range/az looks to look down data. because alos range\n'); 
        fprintf(fid, '       to azimuth resolution is 1:2, always look down by factor of 2-->\n'); 
        fprintf(fid, '       <property name="range looks">3</property>\n'); 
        fprintf(fid, '       <property name="azimuth looks">6</property>\n\n'); 
        fprintf(fid, '	<component name="master">\n'); 
        fprintf(fid, sprintf('       <property name="IMAGEFILE">[''%s'']</property>\n', IMG1)); 
        fprintf(fid, sprintf('       <property name="LEADERFILE">[''%s'']</property>\n', LED1)); 
        fprintf(fid, '<!--            uncomment next line if one SLC is FBS and the other is FBD-->\n'); 
            fprintf(fid, '       <!--<property name="RESAMPLE_FLAG">dual2single</property>-->\n'); 
        fprintf(fid, '       <property name="OUTPUT">master.raw</property>\n'); 
        fprintf(fid, '	</component>\n\n'); 
        fprintf(fid, '	<component name="slave">\n'); 
        fprintf(fid, sprintf('       <property name="IMAGEFILE">[''%s'']</property>\n', IMG2)); 
        fprintf(fid, sprintf('       <property name="LEADERFILE">[''%s'']</property>\n', LED2)); 
            fprintf(fid, '       <!--<property name="RESAMPLE_FLAG">dual2single</property>-->\n'); 
        fprintf(fid, '       <property name="OUTPUT">slave.raw</property>\n'); 
        fprintf(fid, '	</component>\n\n\n'); 
        fprintf(fid, sprintf('       <property name="filter strength">%g</property>\n\n', filter_strength)); 
        fprintf(fid, '	<property name="unwrap">\n'); 
        fprintf(fid, sprintf('       <value>%s</value>\n', unwrap)); 
        fprintf(fid, '	</property>\n'); 
        fprintf(fid, '	<property name="unwrapper name">\n'); 
        fprintf(fid, '       <value>grass</value>\n'); 
        fprintf(fid, '	</property>\n'); 
        fprintf(fid, '	<property name="geocode list">\n'); 
        fprintf(fid, '       <value>["filt_topophase.flat", "los.rdr", "topophase.cor", "topophase.flat", "phsig.cor"]</value>\n'); 
        fprintf(fid, '	</property>\n'); 
        fprintf(fid, '  <property name="geocode bounding box">\n'); 
        fprintf(fid, '       <value>[43.77 43.91 -123.44 -123.30]</value>\n'); 
        fprintf(fid, '  </property>\n\n'); 
        fprintf(fid, '</component>\n'); 
        fprintf(fid, '</insarApp>\n'); 

        fclose(fid); 
        cd ..
    end


