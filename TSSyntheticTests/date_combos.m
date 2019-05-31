% Reformat old .mat files to be in form: 
% params = 
    % dn_all    (n x 1 vector with date numbers of all SAR dates
    % bl_all    (n x 1 vector with baselines of all SAR dates (1st is 0)
    % intcombos (n x 2 matrix, with each element indicating index value of
    %            date in dn_all, and each row in 1 interferogram)

%% Reformat for old matlab scripts
    ccdir = '/data/pmb229/other/clearcuttingTStest/';
    
    % ALOS P222f870
    load([ccdir 'baselines.mat']);
    load([ccdir 'meancor_bl_dates_area2_HH.mat']);
    load('useints.mat'); % from /data/pmb229/isce/p222f870/Lidar_ints/log/

    u = useints.dn;
    d = baselines2.dn_all; 
    b = baselines2.bl_all; 

% ALOS p222f870: make ~daisy chained pairs, 17 ints. 
    ua = [];
    for i = 1:length(u)
        u1 = u(i,1); 
        u2 = u(i,2); 
        i1 = find(u1 == d); 
        i2 = find(u2 == d); 
        ua = [ua; i1 i2];
    end

%     params = struct('dn_all', d, 'bl_all', b, 'intcombos', ua); 
%     save([ccdir 'ALOS_p222f870_params_tscombos.mat'], 'params'); 

    
% ALOS p222f870: make params for old 32 int pairs
    u  = meancor_bl_dates.good_cor_idx; 
    ud = meancor_bl_dates.dateCombos; 
    ua = [];
    for i = 1:length(u)
        x  = ud(u(i),:);
        i1 = find(d == x(1));
        i2 = find(d == x(2));
        ua = [ua; i1 i2];
    end
    
     params = struct('dn_all', d, 'bl_all', b, 'intcombos', ua); 
%     save([ccdir 'ALOS_p222f870_params_tscombos_old.mat'], 'params'); 


% ALOS p222f870: make params for all 210 int pairs
    u  = meancor_bl_dates.dateCombos; 
    ua = [];
    for i = 1:length(u)
        i1 = find(d == u(i,1));
        i2 = find(d == u(i,2));
        ua = [ua; i1 i2];
    end
    
%     params = struct('dn_all', d, 'bl_all', b, 'intcombos', ua); 
%     save([ccdir 'ALOS_p222f870_params_tscombos_all.mat'], 'params'); 
    



% Sentinel data
    load([ccdir 'cctest_sentinel_bl.mat']); 
    load([ccdir 'cctest_sentinel_dc.mat']); 
    d  = bl.dn_all; 
    b  = bl.bl_all; 
    u  = dc.good_cor_idx; 
    ud = dc.dateCombos; 
    ua = [];
    for i = 1:length(u)
        x  = ud(u(i),:);
        i1 = find(d == x(1));
        i2 = find(d == x(2));
        ua = [ua; i1 i2];
    end

    %params = struct('dn_all', d, 'bl_all', b, 'intcombos', ua); 
%     save([ccdir 'Sentinel_params_1.mat'], 'params'); 


    
% Sentinel data - random
    d  = [datenum('01062015', 'ddmmyyyy'):12:datenum('01042019', 'ddmmyyyy')]'; 
    b  = [0; (randn(length(d)-1,1)*70)]; 
    ud = [[1:length(d)-1]' [2:length(d)]']; 
    ud = d(ud); 
    u  = [1:length(ud)]';
    ua = [];
    for i = 1:length(u)
        x  = ud(u(i),:);
        i1 = find(d == x(1));
        i2 = find(d == x(2));
        ua = [ua; i1 i2];
    end

     params2 = struct('dn_all', d, 'bl_all', b, 'intcombos', ua); 
%     save([ccdir 'rand_params_1.mat'], 'params'); 








%% baseline plot
    dn_all = params.dn_all;
    bl_all = params.bl_all;
    ic     = params.intcombos;
    bl     = bl_all(ic); 
    dn     = dn_all(ic); 
    figure; hold on; box on; grid on; 
    for i=1:length(ic)
        plot(dn(i,:), bl(i,:), 'k')
    end
    plot(dn_all, bl_all, '.k', 'markersize', 15);
    datetick; 
    xlabel('date'); ylabel('baseline (m)'); 
    
    
%% baseline plot for andes pres
close; 
figure; hold on; 
yl = [-4200 4000];
subplot(2,1,1); hold on; box on; grid on; 
    dn_all = params.dn_all;
    bl_all = params.bl_all;
    ic     = params.intcombos;
    bl     = bl_all(ic); 
    dn     = dn_all(ic); 
    for i=1:length(ic)
        plot(dn(i,:), bl(i,:), 'k')
    end
    plot(dn_all, bl_all, '.k', 'markersize', 10);
    datetick; 
    ylabel('Baseline (m)'); 
    ylim(yl); 
subplot(2,1,2); hold on; box on; grid on; 
    dn_all = params2.dn_all;
    bl_all = params2.bl_all;
    ic     = params2.intcombos;
    bl     = bl_all(ic); 
    dn     = dn_all(ic); 
    for i=1:length(ic)
        plot(dn(i,:), bl(i,:), 'k')
    end
    plot(dn_all, bl_all, '.k', 'markersize', 10);
    datetick; 
    xlabel('Date'); ylabel('Baseline (m)'); 
    ylim(yl); 
    xlim([datenum('01012015', 'ddmmyyyy') datenum('01012020', 'ddmmyyyy')]);
    
    

    