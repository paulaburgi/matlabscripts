% cctest_vel_changedz.m

% close all
clear 

ccdir = '/data/pmb229/other/clearcuttingTStest/'; 

% load baseline data & int data
%     load([ccdir 'ALOS_p222f870_params_tscombos.mat']); %17 ints
    load([ccdir 'ALOS_p222f870_params_tscombos_old.mat']); %32 ints
%    load([ccdir 'ALOS_p222f870_params_tscombos_all.mat']); %210 ints
    	l = 0.236; sr = 8.39e5; los = 38.7; % ALOS-1
%      load([ccdir 'Radarsat2_francisco_params_1.mat']); 
%         l = 0.055; sr = 9.2e5; los = 34.0; % Radarsat-2
%     load([ccdir 'Sentinel_params_1.mat']); 
%         l = 0.055; sr = 6.9e5; los = 43.0; % Sentinel 
rparams = [l; sr; los]; 

% extract useful info
    dn_all = params.dn_all; 
    bl_all = params.bl_all; 
    gidx_orig   = params.intcombos; 
    dc_orig     = dn_all(gidx_orig); 
    du          = sort(unique(dc_orig(:))); 

% dates of dem error introduction
    ccdates  = [du(1)-15; du+1]; 
    %ccdates  = du(1)-16; 
    nccdates = length(ccdates); 
    ntests   = 1; 
    nw       = 0.0; % nw=1 ~ +/-1m noise (disp = lamda*phase / 4*pi) 


mb_all   = [];
mv_all   = {};
mb_all2  = [];
mv_all2  = {};
% loop through date of DEM introduction, remake G every time
for i = 1:nccdates
    c = ccdates(i); 
    btwint = zeros(length(dc_orig), 1); 
    aftint = ones(length(dc_orig), 1); 
    mb     = [];
    mv     = [];
    mb2    = [];
    mv2    = [];
    % find ints we want to take out
    for j = 1:length(dc_orig)
        idx1 = dc_orig(j,1) <= c;
        idx2 = dc_orig(j,2) >= c;
        if idx1 && idx2
            btwint(j) = 1; 
            aftint(j) = 0; 
        elseif idx1
            aftint(j) = 0; 
        end
    end
    gidx   = gidx_orig(~btwint,:); 
    aftint = aftint(~btwint,:); 
    bl     = bl_all(gidx); 
    bl     = (bl(:,2)-bl(:,1));
    dc     = dn_all(gidx); 
    nvels  = length(du)-1;
    nints  = length(dc); 
    [Gbl, zrdi, blg] = make_cctest_G(dc, du, bl, rparams);
    
    % DEM error in ints
    dz     = 30;                                           % initial topo error
    dz2     = [make_variable_dz(dz, dc, 0.00274); zrdi]; % topo error with tree growth
    
    % define time steps and "real" deformation
    d2y  = 0.00274;        % days to years
    vel  = ones(nvels-1,1)*0; %(diff(du)*d2y)*1; % vel in cm/yr
    def  = Gbl(1:end-length(zrdi),1:end-1)*vel; 
    % create ints
    intsr   = [def; zrdi];
    iz      = blg.*dz;   % addition to phase due to topo error
    iz2     = blg.*dz2;   % addition to phase due to topo error
    intsiz  = intsr;
    intsiz2 = intsr;
    aftint  = [aftint; zrdi]; 
    intsiz(find(aftint == 1))  = intsiz(find(aftint == 1)) + iz(find(aftint == 1)); 
    intsiz2(find(aftint == 1)) = intsiz2(find(aftint == 1)) + iz2(find(aftint == 1)); 
    
    for k = 1:ntests
        n        = [(rand(length(dc), 1)-0.5).*nw; zrdi];
        
        intsizn  = intsiz + n;
        intsizn2 = intsiz2 + n;
        mest     = inv(Gbl'*Gbl)*Gbl'*intsizn; 
        mb       = [mb; mest(end)]; 
        mv       = [mv mest(1:end-1)]; 
        mest2     = inv(Gbl'*Gbl)*Gbl'*intsizn2; 
        mb2       = [mb2; mest2(end)]; 
        mv2       = [mv2 mest2(1:end-1)]; 
    end
    mb_all  = [mb_all; mb]; 
    mv_all  = [mv_all; mv]; 
    mb_all2 = [mb_all2; mb2]; 
    mv_all2 = [mv_all2; mv2]; 
    
    
   
end


%% Plot

% plot full time series with velocity at each timestep 
figure; hold on; box on; 
plot([du(1)-500 du(end)+500], [0 0], 'color', [0.5 0.5 0.5]); 
ylim([-15 15]); 
cmaps = jet(nccdates)*0.8; 
m  = [];
s  = [];
m2 = [];
s2 = [];
for j = 1:nccdates
    mvi  = cell2mat(mv_all(j)); 
    mvi2 = cell2mat(mv_all2(j)); 
    for i = 1:ntests
        plot(du, [0; 0; mvi(:,i)], '-', 'color', cmaps(j,:), 'linewidth', 1); 
         plot(du, [0; 0; mvi2(:,i)], '--', 'color', cmaps(j,:), 'linewidth', 1); 
    end
    m = [m; mean(mean(mvi))]; 
    s = [s; std(mean(mvi))]; 
    m2 = [m2; mean(mean(mvi2))]; 
    s2 = [s2; std(mean(mvi2))]; 
end
datetick; 
xlim([du(1)-100 du(end)+100]);
xlabel('Date'); 
ylabel('Velocity (cm/yr)'); 


% plot mean velocity vs date of DEM error introduction, with standard dev.
figure; hold on; box on;
plot([du(1)-1000 du(end)+1000], [0 0], 'k'); 
errorbar(ccdates, m, s, 'linewidth', 1);
errorbar(ccdates, m2, s2, 'linewidth', 1);
xlabel('Date of DEM introduction'); 
ylabel('Average Velocity'); 
datetick; 
xlim([du(1)-100 du(end)+100]); 

keyboard














