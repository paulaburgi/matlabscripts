%cormask_v1.m

ifol = '/data/pmb229/isce/p222f870/stack/Igrams/';
cd(ifol); 

clear;

% load int info
load('/data/pmb229/isce/p222f870/data/baselines/baselines.mat'); 
load('/data/pmb229/isce/p222f870/data/analysis/meancor_bl_dates_area2_HH.mat'); 
d      = meancor_bl_dates; 
gidx   = d.good_cor_idx; 
dc     = d.dateCombos; 
dc     = dc(gidx,:); 
gidx2  = [0, 3, 5, 8, 11, 14, 16, 19, 20, 22, 23, 24, 26, 27, 28, 30, 31]+1; 
dc     = dc(gidx2, :); 
d1     = [repmat('20', length(dc), 1) datestr(dc(:,1), 'yymmdd')]; 
d2     = [repmat('20', length(dc), 1) datestr(dc(:,2), 'yymmdd')]; 
nints  = length(dc); 

% get nx and ny
fol1 = [d1(1,:) '_' d2(1,:)];
x=importdata([fol1 '/filt_' fol1 '.cor.vrt']);
l1 = x{1}; 
qf = strfind(l1, '"'); 
nx = str2num(l1(qf(1)+1:qf(2)-1)); 
ny = str2num(l1(qf(3)+1:qf(4)-1)); 

% load all cor files
cor = [];
for i = 1:nints
    % load cor
    foli = [d1(i,:) '_' d2(i,:)];
    
    % FIX THIS - ONE SLC IS MISSING
    if exist(foli, 'file'); 
        f    = [foli '/filt_' foli '.cor'];
        fid  = fopen(f, 'r', 'native');
        t    = fread(fid, [nx, ny], 'real*4');
        fclose(fid);
        ii = size(cor, 3); 
        cor(:,:,ii) = t; 
    end
end

lm = 0.4; 
corm = min(cor, [], 3); 
idx = find(corm < lm); 
idx2 = find(corm >= lm); 
corm(idx) = 0; 
%corm(idx2) = 1; 
pcolor(corm(1:2:end,1:2:end)'); shading flat; colormap gray


% reading a msk file is like this: 
%t = fread(fid, [3462, 5290], 'uint8');
pbfol = '/data/pmb229/isce/p222f870/stack/merged/pbints/';
fw    = [pbfol 'meancor.msk']; 
fid   = fopen(fw, 'w'); 
fwrite(fid, fw, 'float32'); 
fclose(fid); 

% kyle says to run snaphu with configuration file to give it a coherence
% threshold (look at online manual)
% run it with snaphu -f int.config
% e.g. https://web.stanford.edu/group/radar/softwareandlinks/sw/snaphu/snaphu.conf.brief

















