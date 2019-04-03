%% Load data, build basic G

close all
%clear 

ccdir = '/data/pmb229/other/clearcuttingTStest/'; 

% load baseline data & int data
%     load([ccdir 'ALOS_p222f870_params_tscombos.mat']); %17 ints
%     load([ccdir 'ALOS_p222f870_params_tscombos_old.mat']); %32 ints
   load([ccdir 'ALOS_p222f870_params_tscombos_all.mat']); %210 ints
    	l = 0.236; sr = 8.5e5; los = 38.7; % ALOS-1
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
    du_orig     = sort(unique(dc_orig(:))); 
    
    %gidx   = [gidx(1:4,:); gidx(6:end,:)]; 

% dates of dem error introduction
    ccdates  = [du_orig(1)-15; du_orig+1]; 
    %ccdates  = du_orig(5)+1; 
    nccdates = length(ccdates); 
    ntests   = 100; 
    nw       = 0.1; % nw=1 ~ +/-1cm noise (disp = lamda*phase / 4*pi) 
    
mb_all   = [];
mv_all   = {};
% loop through date of DEM introduction, remake G every time
for i = 1:nccdates
    c = ccdates(i); 
    btwint = zeros(length(dc_orig), 1); 
    aftint = ones(length(dc_orig), 1); 
    mb     = [];
    mv     = [];
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
    du     = sort(unique(dc_orig(:))); 
    nvels  = length(du)-1;
    nints  = length(dc); 
    [Gbl, zrdi, blg] = make_cctest_G(dc, du, bl, rparams);
    
    % DEM error in ints
    dz     = 30;                                        
    
    % define time steps and "real" deformation
    vel  = ones(nvels-1, 1).*0; 
    def  = Gbl(1:end-length(zrdi),1:end-1)*vel; 
    % create ints
    intsr = [def; zrdi];
    iz     = blg.*dz;  
    intsiz = intsr;
    aftint = [aftint; zrdi]; 
    intsiz(find(aftint == 1)) = intsiz(find(aftint == 1)) + iz(find(aftint == 1)); 
    
    for k = 1:ntests
        n        = [(rand(length(dc), 1)-0.5).*nw; zrdi];

        intsizn  = intsiz + n;
        mest     = inv(Gbl'*Gbl)*Gbl'*intsizn; 
        mb       = [mb; mest(end)]; 
        mv       = [mv mest(1:end-1)]; 
    end
    mb_all = [mb_all; mb]; 
    mv_all = [mv_all; mv]; 
    
    
   
end


%% Plot

figure; hold on; box on; 
ylim([-30 30]); 
cmaps = jet(nccdates); 
m = [];
s = [];
for j = 1:nccdates
    mvi = cell2mat(mv_all(j)); 
    for i = 1:ntests
        plot(du, [0; 0; mvi(:,i)], '*-', 'color', cmaps(j,:)); 
    end
    m = [m; mean(mean(mvi))]; 
    s = [s; std(mean(mvi))]; 
end
datetick; 
xlim([du(1)-100 du(end)+100]);
xlabel('Date'); 
ylabel('Velocity (cm/yr)'); 

figure; hold on; box on;
plot([du(1)-1000 du(end)+1000], [0 0], 'k'); 
errorbar(ccdates, m, s);
xlabel('Date of DEM introduction'); 
ylabel('Average Velocity'); 
datetick; 
xlim([du(1)-100 du(end)+100]); 

keyboard














