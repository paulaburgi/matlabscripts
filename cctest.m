%% Load data, build basic G

close all

% load baseline data
load('baselines.mat');
dn_all = baselines2.dn_all; 
bl_all = baselines2.bl_all; 

% load interferogram combos
load('meancor_bl_dates_area2_HH.mat');
d      = meancor_bl_dates; 
gidx   = d.good_cor_idx; 
dc     = d.dateCombos; 
bl     = d.bl; 

% extra ints to connect with network
exints = [733236 733604; 733558 733604; 733742 733972; 733788 733972];
[a,b]  = dsearchn(dc, exints);
gidx   = [gidx; a]; 

% select which ints to work with
bl     = bl(gidx,:); 
dc     = dc(gidx,:); 
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
nw  = 5; 
n   = (rand(length(dc), 1)-0.5).*nw;

% No DEM error
    % define time steps and "real" deformation
    t     = cumsum(diff(du));
    def   = [1:length(du)-1]'; %zeros(length(du)); 
    intsr = G*def; 
    % add noise to interferograms
    ints  = intsr+n;
    % do inversion
    m     = inv(G'*G)*G'*ints; 

% add baseline term to G matrix
    blg    = (4.*pi.*bl)./(.236.*85e5.*sind(38.7)); 
    Gbl    = [G blg];
    dz     = 30; 
    defz   = [def; dz];
    intsrz = Gbl*defz; 
    intsz  = intsrz+n; 
    mz     = inv(Gbl'*Gbl)*Gbl'*intsz; 

% figure; hold on; 
% plot(t, def, '-*k'); 
% plot(t, m, '-*r'); 
% plot(t, mz(1:end-1), '-*g'); 
% close all

%% Add it DEM error with each time step 

dz  = 30; % topo error
iz  = blg.*dz; 
idz = intsr; 
mbl = [];
nwk = linspace(0,.1,100); 
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
% x label info
    xlab = datestr(str2double(xticklabels).*1e5, 'mm-dd-yyyy'); 
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
    axis([xt(1)-50 xt(end)-60 min(bl_all)-3000 max(bl_all)+3000])
% plot DEM error estimations
subplot(1,3,2)
    hold on; box on; grid on
    plot(du, flipud(mbl_all), '-', 'linewidth', 1, 'color', [.7 .7 .7]); 
    plot(du, flipud(mbl_all(:,1)), '-ko', 'linewidth', 2); 
    plot(du(1), dz, '.', 'markersize', 30, 'color', [0.6 0 0]); 
    set(gca,'xtick',xt); 
    set(gca,'xticklabel',xlab);
    set(gca,'XTickLabelRotation',45);
    xlabel('Date of clear cutting'); 
    ylabel('Estimated DEM error (m)'); 
    title('Temporally Variable DEM Error'); 
    axis([xt(1)-50 xt(end)-60 min(min(mbl_all))-1 dz+1])
    
    ts1   = ts_all(:,1); 
    ts1r  = repmat(ts1, 1, size(ts_all,2));
    dts   = sum((ts1r-ts_all).^2); 
    rmsts = sqrt(dts./length(du)); 
    plot(du, rmsts.*10, '.-', 'markersize', 8, 'color', [.7 0 0]); 
%plot time series for no-noise models
subplot(1,3,3)
    hold on; box on; grid on
    for i=1:size(ts_all,2)
        plot(du, [0; ts_all(:,i)], '-', 'linewidth', 1, 'color', cj(i,:)); 
    end
    cj = jet(length(ts_all)); 
    plot(du, [0; ts_all(:,1)], 'k-', 'linewidth', 3); 
    axis([xt(1)-50 xt(end)-60 min(min(ts_all)) max(max(ts_all))]); 
    set(gca,'xtick',xt); 
    set(gca,'xticklabel',xlab);
    set(gca,'XTickLabelRotation',45);
    xlabel('Date'); 
    ylabel('RMS error (m)'); 
    title('RMS Error of Time Series');






















