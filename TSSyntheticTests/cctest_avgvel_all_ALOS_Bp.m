%% Load data, build basic G

close all
clear 

ccdir = '/data/pmb229/other/clearcuttingTStest/'; 
cd(ccdir); 

% load baseline data & int data
load('ALOS_all_params_d_bl_d150-bl1000.mat');       % _d100-bl800, _d150-bl1000, d100000000-bl100000000
load('ALOS_all_params_gidx_all_d150-bl1000.mat');   % _d100-bl800, _d150-bl1000, d100000000-bl100000000
% load('rand_all_params_d_bl_d150-bl1000.mat'); 
% load('rand_all_params_gidx_all_d150-bl1000.mat'); 

% l = 0.055; sr = 9.48e5; los = 43.0;
l           = 0.236; 
sr          = 8.5e5; 
los         = 38.7;

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
            G      = (diff(dc')'*0.002744*4*pi)/l; 
            blg    = (4*pi*bl)./(l*sr*sind(los));  
            %blg    = (4*pi*bl)./(l*sr*sind(los)); 
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
%         intsiz(aftint == 1) = intsiz(aftint == 1) + iz(aftint == 1); 
        intsiz(aftint == 0) = intsiz(aftint == 0) + iz(aftint == 0);  % 'lidar'
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
%         disp([num2str(round(percentDone*100)) '% done,  time left: ' str]); 
    end
        
   
end

%% dz and vel
close all
figure('units', 'normalized', 'outerposition', [.1 .3 .37 .7]); hold on; box on; linkaxes; 
subplot('position', [0.1 0.53 0.8 0.4]); hold on; box on; 
plot([datenum('2006-01-01') datenum('2013-06-01')], [0 0], 'color', [0.8 0.8 0.8], 'linewidth', 2); 
subplot('position', [0.1 0.1 0.8 0.4]);  hold on; box on; 
%fill([datenum('2006-09-04') datenum('2011-08-30') datenum('2011-08-30') datenum('2006-09-04')], [dz dz 0 0], [0.9 0.9 0.9], 'EdgeColor','none'); 
plot([datenum('2011-08-30') datenum('2006-09-04')], [dz dz], '--', 'color', [0.6 0.8 0.6], 'linewidth', 2); 
dia = []; 
for j = 1:ncombos
    mvi  = cell2mat(mv_all(j))*-100; 
    mbi  = cell2mat(mb_all(j)); 
    di   = cell2mat(md_all(j)); 
    dia  = [dia; di]; 
    subplot('position', [0.1 0.53 0.8 0.4]); hold on; 
    plot(di, mvi, '.-', 'color', [0.5 0.5 0.5], 'linewidth', 1)
    subplot('position', [0.1 0.1 0.8 0.4]); 
    plot(di, mbi, '.-', 'color', [0.0 0.4 0.1])
end
subplot('position', [0.1 0.53 0.8 0.4]); hold on; 
plot(cell2mat(md_all(57)), cell2mat(mv_all(57))*-100, 'k', 'linewidth', 2);
ylabel('Inferred Avg Velocity (cm/yr)'); 
xlim([datenum('2006-09-01') datenum('2011-09-01')]); 
set(gca, 'FontName', 'Arial', 'fontsize', 12, 'xticklabel', '')
ylim([-8 14])

subplot('position', [0.1 0.1 0.8 0.4]); 
plot(cell2mat(md_all(57)), cell2mat(mb_all(57)), 'k', 'linewidth', 2); 
xlabel('Date of DEM Error Removal'); 
ylabel('DEM Error (m)'); 
datetick; 
xlim([datenum('2006-09-01') datenum('2011-09-01')]); 
set(gca, 'FontName', 'Arial', 'fontsize', 12)
ylim([-6 48])


ud  = unique(dia); 
bb  = ones(length(ud), ncombos)*-999;
bb2 = ones(length(ud), ncombos)*-999;
for i = 1:ncombos
    mvi              = cell2mat(mv_all(i))*-100; 
    mbi              = cell2mat(mb_all(i)); 
    di               = cell2mat(md_all(i)); 
    [Intsct, i1, i2] = intersect(ud, di);
    bb(i1, i)        = mvi; 
    bb2(i1, i)       = mbi; 
end

bbm  = []; 
bb2m = [];
for i =1:length(ud)
    bbi   = bb(i,:);
    bb2i  = bb2(i,:); 
    bbm   = [bbm; mean(bbi(find(bbi ~= -999)))];
    bb2m  = [bb2m; mean(bb2i(find(bb2i ~= -999)))];
end
bbm2  = movmean(bbm, 14); 
bb2m2 = movmean(bb2m, 14);
% bbm = conv(bbm, [0 0 1 1 1 1 1 0 0]./5); 
% bb2m = conv(bb2m, [0 0 1 1 1 1 1 0 0]./5); 

subplot('position', [0.1 0.53 0.8 0.4]); hold on; 
% plot(ud, bbm(4:end-5), 'r'); 
plot(ud, bbm2, 'r'); 
subplot('position', [0.1 0.1 0.8 0.4]); 
% plot(ud, bb2m(4:end-5), 'r'); 
plot(ud, bb2m2, 'r'); 


%% for comparing tree growth 
% figure; hold on; box on;
% plot([datenum('2006-01-01') datenum('2013-06-01')], [0 0], 'color', [0.7 0.7 0.7], 'linewidth', 1); 
% cmap = brewermap(ncombos,'Accent'); 
% df = {};
% for j = 1:ncombos
%     mvi = cell2mat(mv_all(j))*-100; 
%     mi  = mean(mvi); 
%     di  = cell2mat(md_all(j)); 
%     %plot(di, mvi, '-', 'color', cmap(j,:), 'markersize', 5)
%     mvi2 = cell2mat(mv_all2(j))*-100; 
%     df   = [df; mvi-mvi2]; 
%     %plot(di, mvi2, '--', 'color', cmap(j,:), 'markersize', 5)
%     plot(di, mvi2-mvi, 'color', cmap(j,:)); 
%     if j == 57
%        dii = di; mvii = mvi; mvii2 = mvi2; 
%     end
% end
% xlabel('Date of DEM introduction'); 
% ylabel('Average Velocity (cm/yr)'); 
% % plot(dii, mvii, 'k-', 'markersize', 5, 'linewidth', 2)
% % plot(dii, mvii2, 'k--', 'markersize', 5, 'linewidth', 2)
% plot(dii, mvii2-mvii, 'k-', 'markersize', 5, 'linewidth', 2)
% datetick; 
% %ylim([-12 22]);
% xlim([datenum('2006-09-01') datenum('2011-09-01')]); 
% 
% 
% figure; hold on; box on;
% plot([datenum('2006-01-01') datenum('2013-06-01')], [0 0], 'color', [0.7 0.7 0.7], 'linewidth', 1); 
% for j = 1:ncombos
%     dfi = cell2mat(df(j)); 
%     di  = cell2mat(md_all(j)); 
%     plot(di, dfi, '-.', 'color', cmap(j,:), 'markersize', 5)
%     if j == 57
%        dii = di; dfi2 = dfi; 
%     end
% end
% xlabel('Date of DEM introduction'); 
% ylabel('Vavg Difference between static dz & variable dz (cm/yr)'); 
%  plot(dii, dfi2, 'k-.', 'markersize', 5, 'linewidth', 2)
% 
% datetick; 
% xlim([datenum('2006-09-01') datenum('2011-09-01')]); 





















