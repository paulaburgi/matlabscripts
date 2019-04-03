% Synthetic time series for DEM error introduced mid-time series. Baselines
% generated randomly, to simulate Sentinel or Nisar-like imaging
% geometries. 

close all
clear 
ccdir = '/data/pmb229/other/clearcuttingTStest/'; 

% Static variables
    % Satellite parameters
    l          = 0.055; 
    sr         = 9.48e5; 
    los        = 43.0; 
    rparams    = [l; sr; los]; 
    % time series duration
    d1         = '20180101';
    d2         = '20190101';
    dn_all     = [datenum(d1, 'yyyymmdd'):12:datenum(d2, 'yyyymmdd')]';
    % time series connectivity
    gidx_orig  = [[1:length(dn_all)-1]' [2:length(dn_all)]'];
    dc_orig    = dn_all(gidx_orig); 
    du         = sort(unique(dc_orig(:))); 
    nvels      = length(du)-1;
    % velocity and DEM error magnitude
    vel1       = ones(nvels-1, 1)*0; 
    dz         = 30; 
    
% for loop paremeters
    % dates of dem error introduction
    ccdates    = du+1; 
    nccdates   = length(ccdates); 
    % number of rand baselines generated per DEM error introduction step
    ntests     = 1000; 
    % misc
    mb_all     = {};
    mv_all     = {};
    timerVal   = tic;

% loop through date of DEM introduction, remake G every time
for i = 1:nccdates
    % for loop parameters
    c      = ccdates(i); 
    btwint = zeros(length(dc_orig), 1); 
    aftint = ones(length(dc_orig), 1); 
    mb     = [];
    mv     = [];
    % need to set 1st and 2nd date = 0 if you are missing first date pair
    if c>du(1) && c<du(2)
        xx = 1;
    else
        xx = 0; 
    end
    
    % find ints with/without DEM error
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
    % take out ints that span date of DEM error introduction
    gidx      = gidx_orig(~btwint,:); 
    aftint    = aftint(~btwint,:); 
    dc        = dn_all(gidx); 
    nints     = length(dc); 
    
    % make G matrix 
    [G, zrdi] = make_cctest_G_randBl(dc, du);
    if xx
        G    = G(1:end-1,2:end); 
        zrdi = []; % only valid if skip one... fix later
        vel  = vel1(2:end); 
    else
        vel = vel1; 
    end
    % define "real" deformation
        def   = G(1:end-length(zrdi),:)*vel; 
        intsr = [def; zrdi];
    
     for k =1:ntests
         % generate random imaging geometry
         bli   = [0; (randn(length(dn_all),1)*70)]; 
         bli   = bli(gidx); 
         bli   = (bli(:,2)-bli(:,1)); 
         % baseline to phase change of DEM error (in meters)
         blgi  = [(2*bli)./(rparams(2).*sind(rparams(3))); zrdi]; 
         Gbl   = [G blgi];
        
        % add noise to ints, if wanted
        nw      = 0.001;    %0.001~0.5cm noise, 0.01~5cm noise
        %n      = [(randn(length(dc), 1)).*nw; zrdi];
        intsiz  = intsr; % + n;
        
        % add baseline-dependent term to applicable ints
        iz      = blgi.*dz;
        aftint  = [aftint; zrdi]; 
        intsiz(aftint == 1) = intsiz(aftint == 1) + iz(aftint == 1); 
        
        % inversion
        mest     = inv(Gbl'*Gbl)*Gbl'*intsiz; 
        
        % ignore anomolously large velocity estimates
        %if std(mest) < 500       
            % add 0 if no int spans first time step
            if  xx
                mest = [0; mest]; 
            end
            mb       = [mb; mest(end)]; 
            mv       = [mv mest(1:end-1)]; 
        %end
    end
    mb_all = [mb_all mb]; 
    mv_all = [mv_all; mv]; 
    
    % approx time estimate for script to run
    if find(i == 1:10:nccdates)
        elapsedTime = toc(timerVal);
        percentDone = (i./nccdates);
        timeLeft    = (elapsedTime/percentDone)-elapsedTime; 
        str         = datestr(seconds(timeLeft),'HH:MM:SS');
        disp([num2str(round(percentDone*100)) '% done,  time left: ' str]); 
    end
   
end


%% Plot

%figure; hold on; box on; 
% ylim([-150 150]); 
% xlim([du(1)-100 du(end)+100]);
% datetick('x', 'keeplimits'); 
% xlabel('Date'); 
% ylabel('Velocity (cm/yr)'); 
% title('no noise'); 
cmaps = jet(nccdates); 
nframes = ntests; 
F = {};
m = [];
s = [];
for j = 1:nccdates
    mvi = cell2mat(mv_all(j)); 
    f = [];
    for i = 1:nframes
%         plot(du, y, '-', 'color', cmaps(j,:), 'linewidth', 1); 
%         drawnow;  
%         F{(nframes*(j-1))+i} = getframe(gcf); 
%         oops
        %plot(du, y, '.', 'color', [0.5 0.5 0.5]);
    end
    si = std(mvi)'; 
    mi = mean(mvi)'; 
    m  = [m; mean(mi)]; 
    s  = [s; mean(si)]; 
%     hist(s2); 
%     oops
end
% writerObj = VideoWriter([ccdir 'TS_video_noiseFree.avi']);
% writerObj.FrameRate = 10;
% open(writerObj);
% for k=1:length(F)
%     frame = cell2mat(F(k));    
%     writeVideo(writerObj, frame);
% end
% close(writerObj);

figure; hold on; box on;
plot([du(1)-1000 du(end)+1000], [0 0], 'k--'); 
errorbar(ccdates, m*100, s*100, 'linewidth', 2);
plot(ccdates, m*100, '.k', 'markersize', 7); 
xlabel('Date of DEM introduction'); 
ylabel('Average Velocity (cm/yr)'); 
datetick; 
xlim([du(1)-100 du(end)+100]); 
%ylim([-2 2]); 
disp(['mean std: ' num2str(mean(s(2:end-2)*100)) 'cm']); 

% keyboard
