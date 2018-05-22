% triplot.m
% makes a 'triangle plot', which plots the correlation of a bunch of
% interferogram pairs. 

clear 
close all

%% Load variables and define directories

% Define directory and data location
  %oregon
    pf_fol  = '/data/pmb229/isce/p222f870/'; 
    cintfol = 'mostcombos/'; pol = 'HH'; 
    % cintfol = 'HVcombos/'; pol = 'HV'; 
    datafol = [pf_fol 'data/']; 
    ds      = [datafol 'baselines/ds_bl-2000m_dl-5e+06m.mat']; 
    bl      = [datafol 'baselines/baselines.mat'];  
  %sumatra
%     pf_fol  = '/data/pmb229/isce/p446f7190_sumatra/'; 
%     cintfol = 'ints_SRTM/';  pol = 'HH'; 
%     datafol = [pf_fol 'data/']; 
%     ds      = [datafol 'baselines/ds_bl-5e+06m_dl-5e+06m.mat']; 
%     bl      = [datafol 'baselines/baselines.mat']; 

% load data
    intfol  = [pf_fol cintfol];
    cd(intfol);
    load(ds); 
    load(bl); baselines = baselines2; 
    intdirs = dir; 
    intdirs = {intdirs.name}; 



    
    
    
%% Define area of interest

% area of interest in pixels (oregon)
    % area 1
    % a=1; 
    % x1 = 420; x2 = 440; 
    % y1 = 360; y2 = 380; 
    % area 2 (shows best coherence values)
    a=2; 
    x1 = 191; x2 = 203;
    y1 = 335; y2 = 358; 
    % area 3
    % a=3; 
    % x1 = 185; x2 = 195;
    % y1 = 85; y2 = 95; 
    
% area of interst in pixels (sumatra)
    % area 1
%     a=1; 
%     x1 = 340; x2 = 360; 
%     y1 = 130; y2 = 150; 
%     a=2; 
%     x1 = 70; x2 = 80;
%     y1 = 350; y2 = 370; 





%% Loop through all cor files 

% define for loop outputs
    forloopidx = 3:length(intdirs); 
    meancor_all=nan(length(forloopidx),1);
    d1=[]; 
    d2=[];
    l=0; 
    
for i=forloopidx
    intdir = cell2mat(intdirs(i)); 
    cd(intdir);
    d1 = [d1; intdir(5:10)];
    d2 = [d2; intdir(12:17)];
    l=l+1; 
    if exist('topophase.cor.geo', 'file')
        % place d1 d2 here if you get an error

        % find size of cor file 
            x=importdata('topophase.cor.geo.vrt');
            l1 = x{1}; 
            qf = strfind(l1, '"'); 
            nx = str2double(l1(qf(1)+1:qf(2)-1)); 
            ny = str2double(l1(qf(3)+1:qf(4)-1)); 

        % import cor file
            filename = 'topophase.cor.geo'; 
            fid         = fopen(filename,'r','native');
            [rmg,~]     = fread(fid,[nx,ny*2],'real*4'); 
                          fclose(fid); 
            mag         = flipud((rmg(1:nx,1:2:ny*2))');
            phsall      = flipud((rmg(1:nx,2:2:ny*2))');
            cor         = phsall(y1:y2, x1:x2); 
                %close
                %figure; hold on; 
                %pcolor(phsall); shading flat; colormap jet; 
                %plot([x1 x1 x2 x2 x1], [y1 y2 y2 y1 y1], 'k-'); 

        % isolate area of interest for correlation
            meancor = mean(mean(cor)); 
            meancor_all(l) = meancor; 
    end
   cd ..
end

     
     
     
     
     
     
%% Find indices of int pairs out of all possible combos, and plug in correlation and baseline values

% define slc dates, find all date and baseline combos
    t1 = dir(datafol); 
    t2 = {t1.name};
    hasdig = regexp(t2, '\d*');
    ld = cellfun('length', hasdig); 
    md = cellfun(@mean, hasdig); 
    slcidx = eq(ld, md); 
    nslc = sum(slcidx); 
    t4 = t2(slcidx); 
    dates = char(t4);
    datesdn = datenum(dates, 'yymmdd'); 
    dcombo = combnk(datesdn,2);
    bcombo = combnk(baselines.bl_all, 2); 

% find dates for each interferogram pair 
    dn1 = datenum(d1, 'yymmdd');
    dn2 = datenum(d2, 'yymmdd');
    datenum_pairs = [dn1 dn2];

% find which pairs were made, out of all possible pairs 
% then, plug in the correlation and baseline values to these indices
    [~, dist]=dsearchn(datenum_pairs, dcombo); 
    dif = find(dist == 0); 
    corall = nan(length(dcombo), 1); 
    corall(dif) = meancor_all; 

    blall  = nan(length(dcombo), 1); 
    bl = bcombo(dif,:); 
    bld = abs(diff(bl'))'; 
    blall(dif) = bld; 

    tblall = nan(length(dcombo), 1); 
    tbl = dcombo(dif,:); 
    tbld = abs(diff(tbl'))'; 
    tblall(dif) = tbld; 

% define "good" date combos
    cor_lim = 0.2; 
    good_cor_idx = find(corall >= cor_lim); 
    
% index of ints on diagonal of triplot
    
    
% define structure to save 
meancor_bl_dates = struct('meancor', corall, 'bl', blall, 'dateCombos', dcombo, ...
                          'good_cor_idx', good_cor_idx);    
    
                      

    
    
    
%% Make square matrices of correlation and baseline for pcolor

% Define indices of the square matrix that will relate the correlation
% value placement
    dn = datesdn; 
    nd=length(dn);
    cids=[];
    
    for i=1:nd
        cids=[cids (nd)*(i-1)+[i:nd-1]];
    end
    
    tri_cor       = nan(nd);     tri_cor(cids)       = corall;
    tri_bl        = nan(nd);     tri_bl(cids)        = blall; 
    tri_datepairs = nan(nd);     tri_datepairs(cids) = tblall; 

    dcids = diff(cids); 
    didx = [1 find(dcids > 1)+1];
    meancor_bl_dates.diagonal_idx = didx; 
    
    
    
    
%% Plot! 

figure('units', 'normalized', 'outerposition', [0 0 1 1]); 
set(gcf, 'PaperOrientation', 'landscape');

% plot 1, colormap
    %ax1 = subplot(1,2,1); hold on; box on; 
    pcolor(dn,dn,tri_cor'); %shading flat;
    cmap = jet; 
    %colormap(ax1, cmap); 
    colormap(cmap)
    h = colorbar;
    ylabel(h, 'Coherence'); 
    caxis([.21 .7]); %caxis([.2 .4]);
    %title(['Oregon ' pol ' Pairs']);
% change axes     
    hyph     = repmat('/', nslc,1); % hyphen
    yr     = repmat('20', nslc, 1); % year
    axlab = [dates(:,3:4) hyph dates(:,5:6) hyph yr dates(:,1:2)]; %label axes mm-dd-yyyy
    axlab = [dates(:,3:4) hyph yr dates(:,1:2)];
    set(gca,'xtick', datesdn); 
    set(gca,'ytick', datesdn); 
    set(gca,'xticklabel',axlab);
    set(gca,'yticklabel',axlab);
    set(gca,'XTickLabelRotation',45);
    set(gca,'YTickLabelRotation',45);
    axis([min(dn) max(dn) min(dn) max(dn)]); 

% plot 2, colormap
    ax2 = subplot(1,2,2); hold on; box on; 
    pcolor(dn, dn, tri_bl'); %shading flat; 
    cmap = brewermap(100,'RdBu'); %cmap=[cmap(10:100,:)]; %RdYlGn
    colormap(ax2, (cmap)); 
    h2 = colorbar;
    ylabel(h2, 'spatial baseline (m)'); 
    caxis([26 3500]); 
% change axes
    set(gca,'xtick', datesdn); 
    set(gca,'ytick', datesdn); 
    set(gca,'xticklabel',axlab);
    set(gca,'yticklabel',axlab);
    set(gca,'XTickLabelRotation',90);
    axis([min(dn) max(dn) min(dn) max(dn)]); 
% link axes
    linkaxes

    
%% Saving and printing 
     
% save data: [number or all possible combos x 1], where data is filled into
% portions of this array that have a value. 
    dosave = 0; % if you want to save, set dosave = 1; 
    if dosave == 1
        save([datafol 'analysis/meancor_bl_dates_area' num2str(a) '_' pol '.mat'], 'meancor_bl_dates');
    end

% save plot
    doprint = 0; 
    if doprint == 1
        print(gcf, [datafol 'analysis/triplot_area' num2str(a) '_' pol '.jpg'], '-djpeg');
    end

    
    
    
     




%% Extra things 
keyboard
% % plot only months
%     datetick('x', 12); 
%     datetick('y', 12);
%     axis([datesdn(1), datesdn(end), datesdn(1), datesdn(end)]); 
%     f=gca; 
%     f.YTick = f.XTick;
%     f.YTickLabel = f.XTickLabel;
% 
% % import int file
    x=importdata('topophase.cor.geo.vrt');
    l1 = x{1}; 
    qf = strfind(l1, '"'); 
    nx = str2num(l1(qf(1)+1:qf(2)-1)); 
    ny = str2num(l1(qf(3)+1:qf(4)-1)); 
    
    im = sqrt(-1); 
    filename2   = 'filt_topophase.flat.geo'; 
    fid         = fopen(filename2, 'r','native');
    [rmg,count] = fread(fid,[nx*2,ny],'real*4');
    status      = fclose(fid);
    real        = flipud((rmg(1:2:nx*2,1:ny))');
    imag        = flipud((rmg(2:2:nx*2,1:ny))');
    mag2         = abs(real+im*imag);
    phs2         = angle(real+im*imag);
    pcolor(phs2); shading flat











