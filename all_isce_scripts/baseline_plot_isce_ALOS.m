% make baseline plot
% must be in data/baselines folder,

clear 

% data folder
  % oregon 
%     datafol = '/data/pmb229/isce/p222f870/data/'; 
%     baselinefol = [datafol 'baselines/']; 
  % sumatra
    datafol = '/data/pmb229/isce/p446f7190_sumatra/data/'; 
    baselinefol = [datafol 'baselines/']; 

    cd(baselinefol)
    warning off

% SET TEMPORAL AND SPATIAL BASELINES
    bl_lim = 800; 
    dt_lim = 100; % 2 years: 730 

    % this is only done to find the baseline infomation, has no affect on which date combos you make
    sz = dir('int_*'); 
    d1_all=[];
    d2_all=[];
    bl_all=[]; 
    dn_all=[];
    
% for loop to process ints to get baselines
    for i=1:length(sz)
        sn = sz(i).name; 
        d1 = sn(5:10);
        d2 = sn(12:17); 
        cd(sn)

        if ~exist('slave.iq.vrt')
            system(sprintf('insarApp.py %s.xml --steps --end=preprocess', sn)); 
        end

        [~,blc] = system('more isce.log | grep perp_baseline_top');
        bl = str2num(blc(29:end));
        bl_all = [bl_all; bl];
        d1_all = [d1_all; d1];
        d2_all = [d2_all; d2];
        dn_all = [dn_all; datenum(d2,'yymmdd')]; 
        cd ..
    end
    d1_all = [d1; d1_all];
    d2_all = [d1; d2_all];
    dn_all = [datenum(d1(1,:), 'yymmdd'); dn_all]; 
    bl_all = [0; bl_all]; 

    baselines = struct('d1_all', d1_all, 'd2_all', d2_all, 'dn_all', dn_all, 'bl_all', bl_all);

    blf = 'baselines.mat';
    blfe = exist(blf); 

    if blfe ~= 2 
       save(blf, 'baselines');
    end

% plot
    close all
    figure; hold on; 
    set(gcf, 'Position', [300, 300, 900, 400]) 

    dcom = combnk(dn_all,2);
    bcom = combnk(bl_all,2);

    diffd = abs(diff(dcom'));
    bc = find(diffd<dt_lim);
    dcom_tb = dcom(bc,:); 
    bcom_tb = bcom(bc,:); 

    diffb = abs(diff(bcom_tb')); 
    dc = find(diffb<bl_lim);
    dcom_td = dcom_tb(dc,:); 
    bcom_td = bcom_tb(dc,:); 

    ds1 = datestr(dcom_td(:,1),'yymmdd');
    ds2 = datestr(dcom_td(:,2),'yymmdd');
    bl1 = bcom_td(:,1); 
    bl2 = bcom_td(:,2); 

    bll = length(bl_all); 

    for i=1:length(dcom_td)
        plot(dcom_td(i,:), bcom_td(i,:), 'linewidth', 1.4, 'color', [.8 .0 .2]); 
    end

    xlab = cell(1, bll); 
    for i=[1:bll]
        plot(dn_all(i),bl_all(i), 'k.', 'markersize', 30); 
        xlab(i)={sprintf('%s-%s-%s',d2_all(i,1:2), d2_all(i,3:4), d2_all(i,5:6))}; 
    end
    DateString = datestr(dcom_tb(:,1),'yymmdd');


    yl = ylim; 

    set(gca,'xtick',dn_all); 
    set(gca,'xticklabel',xlab);
    set(gca,'XTickLabelRotation',45);
    %set(gca,'ytick',[roundn(yl(1),3):1e3:roundn(yl(2),3)]); 
    xlim([min(dn_all(:))-25 max(dn_all(:))+25]); 
    %ylim([min(bl_all(:))-400 max(bl_all(:))+400]);
    %xtick([1:bll])
    %xticklabels(xlab)
    %xtickangle(75); 
    box on; 
    xlabel('Date (yy-mm-dd)'); 
    ylabel('Baseline (m)'); 
    set(gca, 'fontsize', 14); 
    ax = gca;
    ax.YGrid = 'on';
    ax.GridLineStyle = '-';
    legend(sprintf('Temporal Baseline: %g days \nSpatial Baseline:     %g m ', dt_lim, bl_lim), 'location', 'southeast');
    ds = struct('ds1', ds1, 'ds2', ds2, 'bl1', bl1, 'bl2', bl2, 'bl', bl1-bl2);
    ds_name=sprintf('ds_bl-%gm_dl-%gm.mat', bl_lim, dt_lim); 
    ds_exist=exist(ds_name); 

    % [~,frm] = system([['more '] sn(1:end-5) ['/*.slc.rsc | grep FRAME']]); 
    % [~,satt] = system([['more '] sn(1:end-5) ['/*.slc.rsc | grep PLATFORM']]); 
    % f = str2num(frm(20:60));
    % ss = satt(20:60); 
    % TF = isspace(ss); 
    % s =  ss(~TF); 
    %title(sprintf('Platform: %s   Path/Frame: %g/%g', s, p, f)); 

    if ds_exist ~= 2 
        save(ds_name, 'ds'); 
    end


    % print(gcf, 'baseline_example2', '-dpdf', '-painters');




