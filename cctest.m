%% Load data, build basic G

% close all
clear 

ccdir = '/data/pmb229/other/clearcuttingTStest/'; 

% load baseline data & int data
    load([ccdir 'baselines.mat']); l = 0.236; sr = 8.5e5; los = 38.7; 
    load([ccdir 'meancor_bl_dates_area2_HH.mat']); 
%     load([ccdir 'Francisco_rs2_bl.mat']); baselines2 = fs; 
%     load([ccdir 'Francisco_rs2_ints.mat']); meancor_bl_dates = fs2; l = 0.055; sr = 9.2e5; los = 34; 
    dn_all = baselines2.dn_all; 
    bl_all = baselines2.bl_all; 
    d      = meancor_bl_dates; 
    gidx   = d.good_cor_idx; 
    dc     = d.dateCombos; 
    bl     = d.bl; 

% extra ints to connect with network
    exints = [733236 733604; 733558 733604; 733742 733972; 733788 733972];
    [a,b]  = dsearchn(dc, exints);
    gidx   = [gidx; a]; 
    % % select which ints to work with
%     bl     = bl(gidx,:); 
%     dc     = dc(gidx,:); 
    dcv    = dc(:);
    du     = sort(unique(dcv));

% build G matrix
    G  = zeros(length(dc), length(du));
    xa = [];
    ya = [];
    for i = 1:length(G)
        xi     = dc(i,1); 
        yi     = dc(i,2); 
        x = find(du == xi); 
        y = find(du == yi); 
        G(i,x) = -1; 
        G(i,y) =  1; 

        xa = [xa; x];
        ya = [ya; y];
    end
    G = G(:,2:end); 

%% Compare data to inv with and without DEM Error
% noise weight
nw  = 1; 
n   = (rand(length(dc), 1)-0.5).*nw;

% No DEM error
    % define time steps and "real" deformation
    dtint = diff(du); 
    t     = cumsum(dtint);
    def   = [1:length(du)-1]'; %zeros(length(du)); 
    intsr = ya-xa; 
    % add noise to interferograms
    ints  = intsr+n;
    % do inversion
    m     = inv(G'*G)*G'*ints; 

% add baseline term to G matrix
    blg    = (4.*pi.*bl)./(l.*sr.*sind(los)); 
    Gbl    = [G blg];
    dz     = 10; 
    defz   = [def; dz];
    intsrz = Gbl*defz; 
    intsz  = intsrz+n; 
    mz     = inv(Gbl'*Gbl)*Gbl'*intsz; 

figure; hold on; 
plot(t, def, '-*k'); 
plot(t, m, '-*r'); 
plot(t, mz(1:end-1), '-*g'); 
close 

%% Add it DEM error with each time step 

dz  = 30; % topo error
iz  = blg.*dz; 
idz = intsr; 
mbl = [];
nwk = linspace(0,1,100); 
mbl_all = [];
ts_all = [];

% loop through noise weighting
for k = 1:length(nwk)
        mbl = [];
    % loop through date of clear cutting
    for i = 1:length(du)
        n   = (rand(length(dc), 1)-0.5).*nwk(k);
        idz = intsr+n; 
        % loops through data values to add Bp affect
        for j = 1:length(dc)
            % change ints that include dates greater than j
            if xa(j) < i
                idz(j) = intsr(j)+iz(j); 
            end
        end
        
        % do inversion
        mz     = inv(Gbl'*Gbl)*Gbl'*idz; 
        mbl    = [mbl; mz(end)]; 
        
        if k == 1
            ts_all = [ts_all mz(1:end-1)];
        end

    end
    mbl_all = [mbl_all mbl]; % accumulate ints with diff noise
end

% Plot
figure; 
% x label info
%     xlab = datestr(str2double(xticklabels).*1e5, 'mm-dd-yyyy'); 
    xlab = datestr(du, 'mm-dd-yyyy'); 
    xlab = ['07-01-2007'; '01-01-2008'; '07-01-2008'; '01-01-2009'; '07-01-2009'; '01-01-2010'; '07-01-2010'; '01-01-2011'; '07-01-2011'];
    xt   = datenum(xlab, 'mm-dd-yyyy'); 
    xlab = datestr(xt, 'mmmyyyy');
    set(gcf, 'Position', [300, 300, 900, 400])
% plot baseline info
subplot(1,3,1)
    hold on; box on; grid on;
    plot(du, bl_all, 'k.', 'markersize', 15); 
    for i=1:length(intsr)
        bl1i = find(dn_all == dc(i,1)); 
        bl2i = find(dn_all == dc(i,2)); 
        bl1  = bl_all(bl1i); 
        bl2  = bl_all(bl2i); 
        plot(dc(i,:), [bl1 bl2], 'k'); 
    end
    set(gca,'xtick',xt); 
    set(gca,'xticklabel',xlab);
    set(gca,'XTickLabelRotation',45);
    xlabel('Date'); 
    ylabel('Bperp (m)'); 
    title('Baseline Plot');
    %axis([xt(1)-50 xt(end)-60 min(bl_all)-3000 max(bl_all)+3000])
% plot DEM error estimations
subplot(1,3,2);
    hold on; box on; grid on
    yyaxis left
    plot(du, ones(length(du),1).*dz, '--', 'linewidth', 2, 'color', [0.7 0.7 0.7]); 
    plot(du, (mbl_all), '-', 'linewidth', 1, 'color', [.7 .7 .7]); 
    plot(du, (mbl_all(:,1)), '-k', 'linewidth', 2); 
    scatter(du, (mbl_all(:,1)), 30, 1:length(du), 'linewidth', 2); colormap jet; 
    set(gca,'xtick',xt); 
    set(gca,'xticklabel',xlab);
    set(gca,'XTickLabelRotation',45);
    xlabel('Date of clear cutting'); 
    ylabel('Estimated DEM error (m)'); 
    title('Temporally Variable DEM Error'); 
    axis([xt(1)-50 xt(end)-60 min(min(mbl_all))-1 dz+5])
    
    yyaxis right
    x=gca; 
    set(x, 'YColor', [0.7 0 0]); 
    ts1   = def(1:end); 
    ts1r  = repmat(ts1, 1, size(ts_all,2));
    dts   = sum((ts1r-ts_all).^2); 
    rmsts = sqrt(dts./length(du)); 
    plot(du, rmsts, '.-', 'markersize', 8, 'color', [.7 0 0]); 
    ylabel('Time Series RMS error', 'color', [0.7 0 0]);
    ylim([0 max(rmsts).*5]); 
    
%plot time series for no-noise models
subplot(1,3,3)
    hold on; box on; grid on
    plot(du, [0; def], 'k-', 'linewidth', 3); 
    cj = jet(length(ts_all)); 
    for i=1:size(ts_all,2)
        plot(du, [0; ts_all(:,i)], '-', 'linewidth', 1, 'color', cj(i,:)); 
    end
    axis([xt(1)-50 xt(end)-60 min(min(ts_all)) max(max(ts_all))]); 
    set(gca,'xtick',xt); 
    set(gca,'xticklabel',xlab);
    set(gca,'XTickLabelRotation',45);
    xlabel('Date'); 
    ylabel('RMS error (m)'); 
    title('RMS Error of Time Series');











% translate francisco data to format of my data
%     fbl = load('Francisco_radarsat2_bl.txt'); 
%     fs  = struct; 
%     fs.dn_all = datenum(num2str(fbl(:,1)), 'yymmdd'); 
%     fs.bl_all = fbl(:,2); 
%     save('Francisco_rs2_bl', 'fs'); 
% 
%     f    = load('Francisco_radarsat2.txt'); 
%     dn   = [datenum(num2str(f(:,1)), 'yyyymmdd') datenum(num2str(f(:,2)), 'yyyymmdd')]; 
%     dcom = fliplr(combnk(fs.dn_all,2));
%     bcom = fliplr(combnk(fs.bl_all,2));
%     bla  = bcom(:,1)-bcom(:,2); 
%     iall = [];
%     for i =1:length(f)
%         idx  = find(dcom(:,1) == dn(i,1) & dcom(:,2) == dn(i,2)); 
%         iall = [iall; idx]; 
%     end
% 
%     dt = [str2num(datestr(dcom(iall,1), 'yyyymmdd')) str2num(datestr(dcom(iall,2), 'yyyymmdd'))];
% 
%     fs2  = struct; 
%     fs2.dateCombos   = dcom; 
%     fs2.bl           = bla; 
%     fs2.good_cor_idx = iall; 
%     save('Francisco_rs2_ints', 'fs2'); 

