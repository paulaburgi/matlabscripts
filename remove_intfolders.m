% filter interferogram folders by temporal baseline, delete interferogram
% folders greater than a certain baseline

clear

% Interferogram folder
    % intfol    = '/data/pmb229/isce/p222f870/NED_ints/'; 
    % intfol    = '/data/pmb229/isce/p222f870/mostcombos/'; 
    intfol    = '/data/pmb229/isce/p222f870/HVcombos/'; 

% Directory names
    cd(intfol); 
    d = dir('int_*'); 
    intdirs = {d.name};

% Baseline and Temporal Limits
    ddlim = 800;  % date difference limit
    bllim = 2000; % baseline limit

% Empty matrices
    dd_all = []; 
    bl_all = []; 
    ddi    = [];
    dd_gt  = [];
    bli    = [];
    bl_gt  = [];
    dt_gt  = [];

% Loop through directories
for i=1:length(intdirs)
    
    intdir = cell2mat(intdirs(i));
    
    % date difference
    d1     = datenum(intdir(5:10), 'yymmdd');
    d2     = datenum(intdir(12:end), 'yymmdd');
    dd     = d2-d1; 
    dd_all = [dd_all; dd];
    if dd > ddlim 
        ddi   = [ddi; i];
        dd_gt = [dd_gt; dd];
        dt_gt = [dt_gt; str2double([num2str(d1) num2str(d2)])];
    end
    
    % baseline
    [~,blstr] = system(['grep perp_baseline_bottom ' intdir '/isce.log']); 
    bidx      = strfind(blstr, '='); 
    bl        = str2double(blstr(bidx+1:end));
    bl_all    = [bl_all; bl];
    if bl > bllim 
        bli = [bli; i];
        bl_gt = [bl_gt; bl];
    end
    
end

% Define which interferogram folders to remove
    removeidx = unique([bli; ddi]);

% Decide whether or not to remove chosen directories
    m=input('Are you sure you want to remove directories? Y/N: '); 
    if m == 'Y';
        for i=1:length(removeidx)
            intdir = cell2mat(intdirs(i));
            system(['rm -r ' intdir '/']); 
        end
    else
        disp('Directories not removed'); 
    end








% Use 'intersect' to see if confirmed coherent pairs are outside your
% bl/date limits
%     metfile   = '/data/pmb229/isce/p222f870/data/analysis/geotiff_gee/meta_all.csv'; 
%     met       = importdata(metfile); 
%     met       = met.textdata; 
%     gd_all    = [];
%     for i = 2:length(met)
%         m = cell2mat(met(i));
%         gd1    = datenum(m(5:10), 'yymmdd');
%         gd2    = datenum(m(12:17), 'yymmdd');
%         gd     = str2num([num2str(gd1) num2str(gd2)]);
%         gd_all = [gd_all; gd];
%     end
