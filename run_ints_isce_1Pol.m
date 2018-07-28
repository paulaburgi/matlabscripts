% run_ints_isce_1Pol.m 
% must have already run: 
    % unzip_rename_data.m
    % Make_baseline_plot.m
    % make_xml_files_1Pol.m
% note: cannot have any other folders/files in this directory besides int folders
% written by Paula Burgi, Jan 24th, 2018

 
% cd to current int dir that you want to process in 
cd('/data/pmb229/isce/p222f870/mostcombos/'); 
%cd('/data/pmb229/isce/p222f870/HVcombos/'); 
%cd('/data/pmb229/isce/p446f7190_sumatra/ints_SRTM'); 
    %cd('/data/pmb229/isce/p222f870/ALOS2_data/iscecombos/'); 
    %cd('/data/pmb229/isce/p222f870/NED_ints/'); 
clear



% parameters 
    isceapp = 'insarApp.py'; % Alos=insarApp.py, sentinel=topsApp.py 
    endend = '--end=endup';
    %restart = '--start=geocode';
    % dostep = '--dostep=filter'; 
    
    
    
    
    
steps = {'startup'; 'preprocess'; 'verifyDEM'; 'pulsetiming'; 'estimateHeights'; ...
         'mocompath'; 'orbit2sch'; 'updatepreprocinfo'; 'formslc'; 'offsetprf'; ...
         'outliers1'; 'prepareresamps'; 'resamp'; 'resamp_image'; 'mocompbaseline'; ...
         'settopoint1'; 'topo'; 'shadecpx2rg'; 'rgoffset'; 'rg_outliers2'; ...
         'resamp_only'; 'settopoint2'; 'correct'; 'coherence'; 'filter'; 'mask'; ...
         'unwrap'; 'unwrap2stage'; 'geocode'; 'endup'};  
intdirs = dir('int_*'); 
intdirs = {intdirs.name}; 
endidx = find(strcmp(steps, endend(7:end)));
nints = length(intdirs); 
nintsdone = 1; 
[~,nints_done] = system('ls */topophase.cor.geo | wc -l'); 
if isnan(str2double(nints_done)); 
    nints_todo = nints; 
else
    nints_todo = nints - str2double(nints_done); 
end


for i=1:nints 
    intdir = cell2mat(intdirs(i)); 
    cd(intdir)
    
    % print to screen
    fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
    fprintf('\n\nrunning %s, %g out of %g\n\n', intdir, nintsdone ,nints_todo)
    
    % start at correct step
    if exist('restart', 'var')
        system(sprintf('%s %s.xml --steps %s %s', isceapp, intdir, restart, endend))
        nintsdone = sum([nintsdone 1]); 
        
    elseif ~exist('PICKLE', 'dir')
        system(sprintf('%s %s.xml --steps %s', isceapp, intdir, endend))
        nintsdone = sum([nintsdone 1]); 
        
    elseif exist('dostep', 'var'); 
        system(sprintf('%s %s.xml --steps %s %s', isceapp, intdir, dostep))
        nintsdone = sum([nintsdone 1]); 
        
    else
        t = dir('PICKLE');
        [~,didx] = sort([t.datenum]); 
        pname = {t.name}; pname = pname(didx); 
        re=regexp(pname, '\.*'); ie=find(cellfun(@isempty,re));
        pname = pname(ie); 
        laststep = cell2mat(pname(end)); 
        lastidx = find(strcmp(steps, laststep)==1); 
            if lastidx >= endidx 
                % do nothing, either processing is done or has completed up to
                % designated "--end=" step
            else 
        start = ['--start=' steps{lastidx+1}]; 
        system(sprintf('%s %s.xml --steps %s %s', isceapp, intdir, start, endend))
        nintsdone = sum([nintsdone 1]); 
            
        end
    end
    cd ..
    fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
    pause(0.1); 
end














