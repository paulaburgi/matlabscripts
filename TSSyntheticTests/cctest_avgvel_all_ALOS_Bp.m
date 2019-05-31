%% Load data, build basic G

close all
clear 

ccdir = '/data/pmb229/other/clearcuttingTStest/'; 
cd(ccdir); 

% load baseline data & int data
load('ALOS_all_params_d_bl_d150-bl1000.mat');       % _d100-bl800, _d150-bl1000, d100000000-bl100000000
load('ALOS_all_params_gidx_all_d150-bl1000.mat');   % _d100-bl800, _d150-bl1000, d100000000-bl100000000
l = 0.055; sr = 9.48e5; los = 43.0; 
rparams =[l; sr; los];
% note, p222f870 is idx=57

% variables
ncombos  = length(d_bl); 
mb_all   = {};
mv_all   = {};
md_all   = {};
mb_all2  = {};
mv_all2  = {};
timerVal = tic;

% loop through date of DEM introduction, remake G every time
for q = 1:ncombos
    % Bl combo
    di          = cell2mat(d_bl(q));
    dn_all      = di(:,1); 
    bl_all      = di(:,2); 
    gidx_orig   = cell2mat(gidx_all(q)); 
    dc_orig     = dn_all(gidx_orig); 
    du          = sort(unique(dc_orig(:)), 'ascend'); 
    nvels       = length(du)-1;
    vel         = 0; 
    dz          = 30; 
    ccdates     = du+1;
    nccdates    = length(ccdates); 
    mb          = [];
    mv          = [];
    mb2         = [];
    mv2         = [];
    
    for i = 1:nccdates
        c = ccdates(i); 
        btwint = zeros(length(dc_orig), 1); 
        aftint = ones(length(dc_orig), 1); 
        na = []; 
        blii = [];
        % need to set second date = 0 too if you are missing first date pair
        if c>du(1) && c<du(2)
            xx = 1;
        else
            xx = 0; 
        end
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
        %gidx   = gidx_orig; 
        aftint = aftint(~btwint,:); 
        blo     = bl_all(gidx); 
        bl     = (blo(:,2)-blo(:,1));
        dc     = dn_all(gidx); 
        nints  = length(dc); 
        dzi    = make_variable_dz(dz, dc, 0.00274); 
    
        % make G matrix 
            G      = diff(dc')'*0.00274; 
            blg    = (2*bl)./(rparams(2).*sind(rparams(3))); 
            Gbl    = [G blg];

        % define time steps and "real" deformation
        def  = Gbl(1:end,1:end-1)*vel; 
        % create ints
        intsr  = def;
        iz     = blg.*dz;
        iz2    = blg.*dzi; 
        
        %nw     = 0.1;    %0.001~0.5cm noise, 0.01~5cm noise
        %n      = [(randn(length(dc), 1)).*nw];
        intsiz = intsr; % + n;
        intsiz(find(aftint == 1)) = intsiz(find(aftint == 1)) + iz(find(aftint == 1)); 
        intsiz2 = intsr; % + n;
        intsiz2(find(aftint == 1)) = intsiz2(find(aftint == 1)) + iz2(find(aftint == 1)); 

        mest     = inv(Gbl'*Gbl)*Gbl'*intsiz; 
        mb       = [mb; mest(2)]; 
        mv       = [mv; mest(1)]; 
        mest2     = inv(Gbl'*Gbl)*Gbl'*intsiz2; 
        mb2       = [mb2; mest2(2)]; 
        mv2       = [mv2; mest2(1)]; 
    end
    mb_all = [mb_all; mb]; 
    mv_all = [mv_all; mv];
    md_all = [md_all; du]; 
    mb_all2 = [mb_all2; mb2]; 
    mv_all2 = [mv_all2; mv2];
    
    
    if find(q == 1:10:ncombos)
        elapsedTime = toc(timerVal);
        percentDone = (q./ncombos);
        timeLeft    = (elapsedTime/percentDone)-elapsedTime; 
        str         = datestr(seconds(timeLeft),'HH:MM:SS');
        disp([num2str(round(percentDone*100)) '% done,  time left: ' str]); 
    end
        
   
end

% dz and vel
figure('units', 'normalized', 'outerposition', [.3 .3 .3 .7]); hold on; box on;
subplot(2,1,1); hold on; box on; 
plot([datenum('2006-01-01') datenum('2013-06-01')], [0 0], 'color', [0.7 0.7 0.7], 'linewidth', 1); 
subplot(2,1,2); hold on; box on; 
% plot([datenum('2006-01-01') datenum('2013-06-01')], [0 0], 'color', [0.7 0.7 0.7], 'linewidth', 1); 
fill([datenum('2006-09-04') datenum('2011-08-30') datenum('2011-08-30') datenum('2006-09-04')], [dz dz 0 0], [0.9 0.9 0.9], 'EdgeColor','none'); 
%cmap = brewermap(ncombos,'Accent'); 
cmap = jet(ncombos)*0.8; 
df = {};
for j = 1:ncombos
    mvi  = cell2mat(mv_all(j))*100; 
    mbi  = cell2mat(mb_all(j)); 
    di   = cell2mat(md_all(j)); 
    subplot(2,1,1); hold on; 
    plot(di, mvi, '-', 'color', cmap(j,:), 'linewidth', 1)
    subplot(2,1,2); hold on; 
    plot(di, mbi, '.-', 'color', cmap(j,:), 'markersize', 5)
end
subplot(2,1,1); 
xlabel('Date of DEM introduction'); 
ylabel('Average Velocity (cm/yr)'); 
datetick; 
xlim([datenum('2006-09-01') datenum('2011-09-01')]); 
%ylim([-12 22]); 
title('150 days, 1000m'); 
subplot(2,1,2); 
xlabel('Date of DEM introduction'); 
ylabel('DEM error (m)'); 
datetick; 
xlim([datenum('2006-09-01') datenum('2011-09-01')]); 




% for comparing tree growth 
figure; hold on; box on;
plot([datenum('2006-01-01') datenum('2013-06-01')], [0 0], 'color', [0.7 0.7 0.7], 'linewidth', 1); 
cmap = brewermap(ncombos,'Accent'); 
df = {};
for j = 1:ncombos
    mvi = cell2mat(mv_all(j))*100; 
    mi  = mean(mvi); 
    di  = cell2mat(md_all(j)); 
    %plot(di, mvi, '-', 'color', cmap(j,:), 'markersize', 5)
    mvi2 = cell2mat(mv_all2(j))*100; 
    df   = [df; mvi-mvi2]; 
    %plot(di, mvi2, '--', 'color', cmap(j,:), 'markersize', 5)
    plot(di, mvi2-mvi, 'color', cmap(j,:)); 
    if j == 57
       dii = di; mvii = mvi; mvii2 = mvi2; 
    end
end
xlabel('Date of DEM introduction'); 
ylabel('Average Velocity (cm/yr)'); 
% plot(dii, mvii, 'k-', 'markersize', 5, 'linewidth', 2)
% plot(dii, mvii2, 'k--', 'markersize', 5, 'linewidth', 2)
plot(dii, mvii2-mvii, 'k-', 'markersize', 5, 'linewidth', 2)
datetick; 
%ylim([-12 22]);
xlim([datenum('2006-09-01') datenum('2011-09-01')]); 


figure; hold on; box on;
plot([datenum('2006-01-01') datenum('2013-06-01')], [0 0], 'color', [0.7 0.7 0.7], 'linewidth', 1); 
for j = 1:ncombos
    dfi = cell2mat(df(j)); 
    di  = cell2mat(md_all(j)); 
    plot(di, dfi, '-.', 'color', cmap(j,:), 'markersize', 5)
    if j == 57
       dii = di; dfi2 = dfi; 
    end
end
xlabel('Date of DEM introduction'); 
ylabel('Vavg Difference between static dz & variable dz (cm/yr)'); 
 plot(dii, dfi2, 'k-.', 'markersize', 5, 'linewidth', 2)

datetick; 
xlim([datenum('2006-09-01') datenum('2011-09-01')]); 





















