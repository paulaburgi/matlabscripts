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
    xmli = [ifol '.xml']; 
    cd(ifol); 
    
    if ~exist('unw2pi.flat', 'file'); 
        % add 2pi increments to non-filtered int
        system(['imageMath.py -e=''round((a_1-arg(b))/(2*PI))'' -o unw2pi.flat ' ...
                '--a=filt_topophase_' fs2 '.unw --b=filt_topophase_' fs2 '.flat']); 
        system(['imageMath.py -e=''a+arg(b)'' -o topophase_' fs2 '.unw ' ...
                '--a=unw2pi.flat --b=topophase.flat']); 
    end
    
    if ~exist(['topophase_' fs2 '.unw.geo'], 'file'); 
        % add topophase_07.unw to geocode file
        system(['sed -i -e ''/geocode list/{n;d}'' ' xmli]);  % removes line after geocode list
        system(['sed -i ''s|<property name="geocode list">|' ... 
                 '<property name="geocode list">\n       <value>' ...
                 '["topophase_' fs2 '.unw"]</value>|g'' ' xmli]); 
        % run isce
        system(['insarApp.py ' xmli ' --steps --start=geocode']); 
    end 
    
    cd ..
end





















