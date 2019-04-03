%% Load data, build basic G

close all
clear 

ccdir = '/data/pmb229/other/clearcuttingTStest/'; 

% load baseline data & int data

dset = 'D1'; % D1=OR ALOS, D2=Fran's RS2, D3=SF Sentinel 

if strcmp(dset, 'D1')
    load([ccdir 'baselines.mat']); 
    load([ccdir 'meancor_bl_dates_area2_HH.mat']); 
    l = 0.236; sr = 8.5e5; los = 38.7; 
elseif strcmp(dset, 'D2')
    load([ccdir 'Francisco_rs2_bl.mat']); baselines2 = fs; 
    load([ccdir 'Francisco_rs2_ints.mat']); meancor_bl_dates=fs2; 
    l=0.055; sr=9.2e5; los=34; 
elseif strcmp(dset, 'D3')
    load([ccdir 'cctest_sentinel_bl.mat']); baselines2 = bl; 
    load([ccdir 'cctest_sentinel_dc.mat']); meancor_bl_dates=dc; 
    l=0.055; sr=6.9e5; lo =43; 
end
dn_all = baselines2.dn_all; 
    bl_all = baselines2.bl_all; 
    d      = meancor_bl_dates; 
    gidx   = d.good_cor_idx; 
    dc     = d.dateCombos; 
    bl     = abs(d.bl); % have to make bl positive for inversion to work
    dta    = diff(dn_all); 

    if strcmp(dset, 'D1')
    % extra ints to connect with network
        exints = [733236 733604; 733742 733972; 733788 733972];
        [a,b]  = dsearchn(dc, exints);
%         gidx   = [gidx; a]; 
    end
    
% % select which ints to work with
    bl     = bl(gidx,:); 
    dc     = dc(gidx,:); 
    dcv    = dc(:);
    du     = sort(unique(dcv));

% build G matrix
    G  = zeros(length(dc), length(du)-1);
    xa = [];
    ya = [];
    for i = 1:length(G)
        xi     = dc(i,1); 
        yi     = dc(i,2); 
        x = find(du == xi); 
        y = find(du == yi); 
        
        dtj   = dta(x:y-1);
        for j = 1:length(dtj)
            G(i,x+j-1) = dtj(j); 
        end

        xa = [xa; x];
        ya = [ya; y];
    end
    Gb  = G; 
    d2y = 0.00274; % convert days to years
    G   = G.*d2y; 

%% Compare data to inv with and without DEM Error
% noise weight
nw  = 0.0; 
n   = (rand(length(dc), 1)-0.5).*nw;

% No DEM error
    % define time steps and "real" deformation
    dtint = diff(du); 
    t     = cumsum(dtint);
    def   = [0:length(du)-1]'; 
    def   = zeros(length(du),1); 
    vel   = diff(def)./dtint; 
    
    
    intsr = def(ya)-def(xa); 
    % add noise to interferograms
    ints  = intsr+n;
    % do inversion
    %m     = inv(G'*G)*G'*ints; 
    % SVD inversion
    r = rank(G); 
    [U,S,V] = svd(G); 
    Gsvd    = V(:,1:r)*inv(S(1:r,1:r))*U(:,1:r)'; 
    msvd    = Gsvd*ints; 
    

% add baseline term to G matrix
    blg    = (4.*pi.*bl)./(l.*sr.*sind(los)); 
    Gbl    = [G abs(blg)]; 
    dz     = 30; 
    velz   = [vel; dz];
    intsrz = intsr+abs(blg.*dz); 
    intsz  = intsrz+n; 
    mz     = inv(Gbl'*Gbl)*Gbl'*intsz; 
    % SVD inversion
    r2      = rank(Gbl); 
    [U,S,V] = svd(Gbl); 
    Gsvd    = V(:,1:r2)*inv(S(1:r2,1:r2))*U(:,1:r2)'; 
    mzsvd   = Gsvd*intsz; 
    
    % vel to def
    %defm     = t.*m; 
    %defmz    = t.*mz(1:end-1); 
    defmsvd  = cumsum(dtint.*msvd); 
    defmzsvd = cumsum(dtint.*mzsvd(1:end-1)); 
    
figure; hold on; box on;
plot(t, def(2:end), '-*k'); 
%plot(t, defm, '-*r'); 
plot(t,defmsvd, '-*b'); 
%plot(t, defmz, '-*g'); 
plot(t, defmzsvd, '-*y'); 
close all

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
        %mz     = inv(Gbl'*Gbl)*Gbl'*idz; 
        mz     = Gsvd*idz; 
        mbl    = [mbl; mz(end)]; 
        
        if k == 1
            ts_all = [ts_all cumsum(dtint.*mz(1:end-1))];
        end

    end
    mbl_all = [mbl_all mbl]; % accumulate ints with diff noise
end





% Plot
% x label info
    dall   = du(1)-100:du(end)+150; 
    [dvec] = datevec(dall); 
    idx1   = find(dvec(:,2) == 1 & dvec(:,3) == 1); 
    idx7   = find(dvec(:,2) == 7 & dvec(:,3) == 1); 
    xt     = sort([dall(idx1) dall(idx7)]); 

set(gcf, 'Position', [300, 300, 900, 500])  
set(gcf, 'PaperOrientation', 'landscape');

% plot baseline info
    pos1 = [0.08 0.18 0.25 0.75];
    subplot('Position',pos1)
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
    set(gca,'XTickLabelRotation',45);
    datetick('x','mm-yyyy', 'keepticks');
    xlabel('Date'); 
    ylabel('Bperp (m)'); 
    title('Baseline Plot');
    dyl = diff(ylim); 
    axis([xt(1)-50 xt(end)-60 min(bl_all)-(dyl/2) max(bl_all)+(dyl/2)])
    
% plot time series for no-noise models 
    pos1 = [0.4 0.381 0.25 0.55];
    subplot('Position',pos1)
    hold on; box on; grid on
    plot(du, [def], 'k-', 'linewidth', 4);
    cj = parula(length(ts_all)); 
    for i=1:size(ts_all,2)
        plot(du, [0;ts_all(:,i)], '-', 'linewidth', 1, ...
            'color', cj(i,:)); 
    end
    axis([xt(1)-50 xt(end)-60 min(min(ts_all)) max(max(ts_all))]); 
    ylabel('deformation (m)');
    title('Time Series & RMS error');
    set(gca,'xtick',xt); 
    set(gca,'xticklabel',{[]}) 
    xl = get(gca, 'xlim'); 
 % plot RMS error   
    pos1 = [0.4 0.18 0.25 0.2];
    subplot('Position',pos1)
    ts1   = def(2:end); 
    ts1r  = repmat(ts1, 1, size(ts_all,2));
    dts   = sum((ts1r-ts_all).^2); 
    rmsts = sqrt(dts./length(du)); 
    plot(du, rmsts, '^-', 'markersize', 4, 'color', [.7 0 0], ...
        'MarkerFaceColor',[.7 0 0]); 
    ylabel('RMS error') %, 'color', [0.7 0 0]);
    set(gca,'xtick',xt); 
    set(gca,'XTickLabelRotation',45);
    datetick('x','mm-yyyy', 'keepticks');
    xlabel('Date'); 
    grid on; 
    xlim([xl]); 

    
%plot DEM error estimations
    pos1 = [0.72 0.18 0.25 0.75];
    subplot('Position',pos1)
    hold on; box on; grid on
    plot(du, ones(length(du),1).*dz, '--', 'linewidth', 2, ... 
        'color', [0.7 0.7 0.7]); 
    plot(du, (mbl_all), '-', 'linewidth', 1, 'color', [.7 .7 .7]); 
    plot(du, (mbl_all(:,1)), '-k', 'linewidth', 2); 
    scatter(du, (mbl_all(:,1)), 30, 1:length(du), 'filled', ...
        'linewidth', 2); 
    colormap parula; 
    set(gca,'xtick',xt); 
    set(gca,'XTickLabelRotation',45);
    datetick('x','mm-yyyy', 'keepticks');
    xlabel('Date of clear cutting'); 
    ylabel('Estimated DEM error (m)'); 
    title('Temporally Variable DEM Error'); 
    axis([xt(1)-50 xt(end)-60 min(min(mbl_all))-1 dz+5])
    
    
    











% translate francisco data to format of my data
%     fbl  = load('Francisco_radarsat2_bl.txt'); 
%     f    = load('Francisco_radarsat2.txt'); 
%     fs   = struct; 
%     dn_all = datenum(num2str(fbl(:,1)), 'yymmdd'); 
%     [dn_all,si] = sort(dn_all); 
%     bl_all = fbl(si,2); 
%     dn   = [datenum(num2str(f(:,1)), 'yyyymmdd') ... 
%            datenum(num2str(f(:,2)), 'yyyymmdd')]; 
%     a    = sort(unique(dn(:)));
%     idx  = ismember(dn_all, a);  
%     fs.dn_all = dn_all(idx); 
%     fs.bl_all = bl_all(idx); 
%     save('Francisco_rs2_bl', 'fs'); 
%     
% 
%     dcom = (combnk(fs.dn_all,2));
%     bcom = (combnk(fs.bl_all,2));
%     bla  = bcom(:,2)-bcom(:,1); 
%     iall = [];
%     for i =1:length(f)
%         idx  = find(dcom(:,1) == dn(i,2) & dcom(:,2) == dn(i,1)); 
%         iall = [iall; idx]; 
%     end
% 
%     dt = [str2num(datestr(dcom(iall,1), 'yyyymmdd')) ... 
%            str2num(datestr(dcom(iall,2), 'yyyymmdd'))];
% 
%     fs2  = struct; 
%     fs2.dateCombos   = dcom; 
%     fs2.bl           = bla; 
%     fs2.good_cor_idx = iall; 
%     save('Francisco_rs2_ints', 'fs2'); 

