% make_xml_files_1Pol.m 
% start in directory that you want all your int_d1_d2 folders to be in 
% makes subfolders within each int_d1_d2 folders, with polarizations
% writes .xml file for each subfolder 
% links data to each subfolder 

% Paula Burgi 
% Written Jan 24, 2018


clear
close all

load('/data/pmb229/isce/p222f870/data/baselines/baselines.mat'); 
load('/data/pmb229/isce/p222f870/data/analysis/meancor_bl_dates_area2_HH.mat'); 

% data folder
    pffol       = '/data/pmb229/isce/p222f870/'; 
    datafol     = [pffol 'data/']; 
    intfol      = [pffol 'mostcombos/']; 

cd(intfol); 

d      = meancor_bl_dates; 
gidx   = d.good_cor_idx; 
dc     = d.dateCombos; 
dc     = dc(gidx,:); 
d1     = datestr(dc(:,1), 'yymmdd'); 
d2     = datestr(dc(:,2), 'yymmdd'); 
nints  = length(dc); 

%filter strength 
fs   = '0.7'; 
fs2  = '07'; 
ofs  = '0.3'; % old filter strength 
ofs2 = '03'; 

for i = 1:nints
    di1   = d1(i, :); 
    di2   = d2(i, :); 
    ifol  = ['int_' di1 '_' di2]; 
    cd(ifol); 
    
    [~, fl] = system(['more ' ifol '.xml | grep filter']);
    
    % copy old files to new name with old filter stregth
    if contains(fl, ofs)
        system(['cp filt_topophase.flat filt_topophase_' ofs2 '.flat']); 
        system(['cp filt_topophase.flat.xml filt_topophase_' ofs2 '.flat.xml']); 
        system(['cp filt_topophase.flat.vrt filt_topophase_' ofs2 '.flat.vrt']);
        system(['cp filt_topophase.unw filt_topophase_' ofs2 '.unw']); 
        system(['cp filt_topophase.unw.vrt filt_topophase_' ofs2 '.unw.vrt']); 
        system(['cp filt_topophase.unw.xml filt_topophase_' ofs2 '.unw.xml']); 
    end
    
    
    if ~contains(fl, fs)      
        
        % replace filter strength 
        xmli = [ifol '.xml']; 
        system(['sed -i ''s|>0.*<|>' fs '<|g'' ' xmli]); 

        % print to screen
        fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
        fprintf('running %s, %g out of %g\n\n', ifol, i ,nints)
        
        % run isce
        system(['insarApp.py ' xmli ' --steps --start=filter --end=unwrap2stage']); 
        % copy new files to new filter strength in name
        system(['cp filt_topophase.flat filt_topophase_' fs2 '.flat']); 
        system(['cp filt_topophase.flat.xml filt_topophase_' fs2 '.flat.xml']); 
        system(['cp filt_topophase.flat.vrt filt_topophase_' fs2 '.flat.vrt']);
        system(['cp filt_topophase.unw filt_topophase_' fs2 '.unw']); 
        system(['cp filt_topophase.unw.vrt filt_topophase_' fs2 '.unw.vrt']); 
        system(['cp filt_topophase.unw.xml filt_topophase_' fs2 '.unw.xml']); 

    end
    
    pause(1)
    
    cd ..
end



