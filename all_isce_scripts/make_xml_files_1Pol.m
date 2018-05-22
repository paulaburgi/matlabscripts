% make_xml_files_1Pol.m 
% start in directory that you want all your int_d1_d2 folders to be in 
% makes subfolders within each int_d1_d2 folders, with polarizations
% writes .xml file for each subfolder 
% links data to each subfolder 

% Paula Burgi 
% Written Jan 24, 2018


clear
close all

% data folder
  % oregon 
    pffol       = '/data/pmb229/isce/p222f870/'; 
    datafol     = [pffol 'data/']; 
    baselinefol = [datafol 'baselines/']; 
    %intfol      = [pffol 'mostcombos/']; 
    %intfol      = [pffol 'HVcombos/']; 
        intfol   = [pffod 'iscecombos']; 
  % sumatra
%     pffol = '/data/pmb229/isce/p446f7190_sumatra/'; 
%     datafol     = [pffol 'data/']; 
%     baselinefol = [datafol 'baselines/']; 
%     intfol = [pffol 'ints_SRTM/']; 

    cd(intfol); 

% matlab file containing the desired date/baseline combinations. 
% e.g. /data/pmb229/roipac/p222f870/data/baselines/
    %datecombofile = 'ds_bl-5e+06m_dl-5e+06m_HV.mat'; %oregon
    %datecombofile = 'ds_bl-5e+06m_dl-500m.mat'; %oregon
    datecombofile = 'ds_bl-2000m_dl-5e+06m.mat';
    ds = [baselinefol datecombofile]; 
    load(ds); 
    ds1 = ds.ds1; 
    ds2 = ds.ds2; 
    nslc = length(ds1); 
    
% parameters
    filter_strength = 0.3;  % 0.0 - 1.0
    unwrap = 'False'; % 'True' or 'False'
     %geocodebox = '[43.77 43.91 -123.44 -123.30]'; % p222f870, orgeon 
     geocodebox = '[43.58 44.38 -123.76 -122.96]'; % p222f870, orgeon
     %geocodebox = '[0.44 0.63 100.37 100.79]'; % p446 f 7190, sumatra

% for loop to make xml files and folders
    for i = 1:nslc
        clear s1 
        clear s2
        s1 = num2str(ds1(i,:)); 
        s2 = num2str(ds2(i,:)); 
        intfol = sprintf('int_%s_%s', s1, s2); 

        % make int folder 
        if exist(intfol) ~= 7
            system(sprintf('mkdir %s', intfol)); 
        end
        cd(intfol)

        % find out if slc are FBD or FBS 
        datafol_s1 = sprintf('%s/%s', datafol, s1); 
        datafol_s2 = sprintf('%s/%s', datafol, s2); 
        s11 = dir(datafol_s1); 
        s21 = dir(datafol_s2); 
        s12 = {s11.name}; 
        s22 = {s21.name}; 
        s13 = cell2mat(s12); 
        s23 = cell2mat(s22); 
            if strfind(s13,'HH') & strfind(s13,'HV')
                s1_fb = 'FBD'; 
            else
                s1_fb = 'FBS'; 
            end
            if strfind(s23,'HH') & strfind(s23,'HV')
                s2_fb = 'FBD'; 
            else
                s2_fb = 'FBS'; 
            end

        % make HHHH/HHHV/HVHV/HVHH folders 
        if strcmp(s1_fb, 'FBS') & strcmp(s2_fb, 'FBS')
            IMG1_FB = 0; 
            IMG2_FB = 0;
        end
        if strcmp(s1_fb, 'FBS') & strcmp(s2_fb, 'FBD')
            IMG1_FB = 0; 
            IMG2_FB = 1; 
        end
        if strcmp(s1_fb, 'FBD') & strcmp(s2_fb, 'FBS') 
            IMG1_FB = 1; 
            IMG2_FB = 0; 
        end
        if strcmp(s1_fb, 'FBD') & strcmp(s2_fb, 'FBD')
            IMG1_FB = 1; 
            IMG2_FB = 1; 
        end


        % define which IMG files to use 
        IMG1 = cell2mat(s12(contains(s12, 'IMG-HH'))); 
        IMG2 = cell2mat(s22(contains(s22, 'IMG-HH'))); 
        LED1 = cell2mat(s12(contains(s12, 'LED'))); 
        LED2 = cell2mat(s22(contains(s22, 'LED'))); 


        %link data files
        system(sprintf('ln -sf %s/%s/%s %s', datafol, s1, IMG1, IMG1)); 
        system(sprintf('ln -sf %s/%s/%s %s', datafol, s2, IMG2, IMG2)); 
        system(sprintf('ln -sf %s/%s/%s %s', datafol, s1, LED1, LED1)); 
        system(sprintf('ln -sf %s/%s/%s %s', datafol, s2, LED2, LED2));

        % write file 
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
        fprintf(fid, '       <!--<property name="posting">30</property>-->\n'); 
        fprintf(fid, '       <property name="range looks">3</property>\n'); 
        fprintf(fid, '       <property name="azimuth looks">6</property>\n\n'); 
        fprintf(fid, '	<component name="master">\n'); 
        fprintf(fid, sprintf('       <property name="IMAGEFILE">[''%s'']</property>\n', IMG1)); 
        fprintf(fid, sprintf('       <property name="LEADERFILE">[''%s'']</property>\n', LED1)); 
        fprintf(fid, '<!--            uncomment next line if one SLC is FBS and the other is FBD-->\n'); 

        if IMG1_FB == 0
            fprintf(fid, '       <!--<property name="RESAMPLE_FLAG">dual2single</property>-->\n'); 
        else
            fprintf(fid, '       <property name="RESAMPLE_FLAG">dual2single</property>\n'); 
        end

        fprintf(fid, '       <property name="OUTPUT">master.raw</property>\n'); 
        fprintf(fid, '	</component>\n\n'); 
        fprintf(fid, '	<component name="slave">\n'); 
        fprintf(fid, sprintf('       <property name="IMAGEFILE">[''%s'']</property>\n', IMG2)); 
        fprintf(fid, sprintf('       <property name="LEADERFILE">[''%s'']</property>\n', LED2)); 

        if IMG2_FB == 0
            fprintf(fid, '       <!--<property name="RESAMPLE_FLAG">dual2single</property>-->\n'); 
        else
            fprintf(fid, '       <property name="RESAMPLE_FLAG">dual2single</property>\n'); 
        end

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
        fprintf(fid, sprintf('       <value>%s</value>\n', geocodebox)); 
        fprintf(fid, '  </property>\n\n'); 
        fprintf(fid, '</component>\n'); 
        fprintf(fid, '</insarApp>\n'); 

        fclose(fid); 
        cd ..
    end





