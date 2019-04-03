%% Load data, build basic G

%close all
%clear 

ccdir = '/data/pmb229/other/clearcuttingTStest/'; 

% load baseline data & int data
load('ALOS_all_params_d_bl_d150-bl1000.mat');       % _d100-bl800, _d150-bl1000
load('ALOS_all_params_gidx_all_d150-bl1000.mat');   % _d100-bl800, _d150-bl1000
l = 0.055; sr = 9.48e5; los = 43.0; 
rparams =[l; sr; los];
% note, p222f870 is idx=57

% variables
ncombos  = length(d_bl); 
mb_all   = {};
mv_all   = {};
md_all   = {};
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
    vel1        = ones(nvels-1, 1)*0; 
    dz          = 30; 
    ccdates     = du+1;
    nccdates    = length(ccdates); 
    mb          = [];
    mv          = [];
    
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
    
        % make G matrix 
        [Gbl, zrdi, blg] = make_cctest_G(dc, du, bl, rparams);
%         if xx
%             Gbl = Gbl(1:end-1,2:end); zrdi = []; vel = vel1(2:end); 
%         else
            vel = vel1; 
%         end

        % define time steps and "real" deformation
        def  = Gbl(1:end-length(zrdi),1:end-1)*vel; 
        % create ints
        intsr = [def; zrdi];
        iz     = blg.*dz;
        
        %nw     = 0.1;    %0.001~0.5cm noise, 0.01~5cm noise
        n      = [(randn(length(dc), 1)).*nw; zrdi];
        intsiz = intsr; % + n;
        aftint = [aftint; zrdi]; 
        intsiz(find(aftint == 1)) = intsiz(find(aftint == 1)) + iz(find(aftint == 1)); 

        Gbl2     = Gbl; %[Gbl(:,1:end-1)]; 
        mest     = inv(Gbl2'*Gbl2)*Gbl2'*intsiz; 
        mb       = [mb; mest(end)]; 
        mv       = [mv mest(1:end-1)]; %-1
    end
    mb_all = [mb_all mb]; 
    mv_all = [mv_all; mv];
    md_all = [md_all; du]; 
    
    
    if find(q == 1:10:ncombos)
        elapsedTime = toc(timerVal);
        percentDone = (q./ncombos);
        timeLeft    = (elapsedTime/percentDone)-elapsedTime; 
        str         = datestr(seconds(timeLeft),'HH:MM:SS');
        disp([num2str(round(percentDone*100)) '% done,  time left: ' str]); 
    end
        
   
end

figure; hold on; box on;
plot([datenum('2006-01-01') datenum('2013-06-01')], [0 0], 'color', [0.7 0.7 0.7], 'linewidth', 2); 
cmap = brewermap(ncombos,'Accent'); 
s_all = [];
for j = 1:ncombos
    mvi = cell2mat(mv_all(j)); 
    mi  = mean(mvi); 
    si  = std(mvi); 
    di  = cell2mat(md_all(j)); 
    %plot(di, mi*100, '.-', 'color', cmap(j,:), 'markersize', 5)
    errorbar(di, mi*100, si*100); 
    s_all = [s_all; si]; 
end
xlabel('Date of DEM introduction'); 
ylabel('Average Velocity (cm/yr)'); 
datetick; 
xlim([datenum('2006-09-01') datenum('2011-09-01')]); 























