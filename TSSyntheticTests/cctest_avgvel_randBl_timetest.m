% Synthetic time series for DEM error introduced mid-time series. Baselines
% generated randomly, to simulate Sentinel or Nisar-like imaging
% geometries. 


%  keyboard; % DON'T DELETE???
%  keyboard; 

%close all
clear 
ccdir = '/data/pmb229/other/clearcuttingTStest/'; 
cd(ccdir); 




%% Set variables 

% Static variables
    % Satellite parameters
    l          = 0.055; 
    sr         = 9.48e5; 
    los        = 43.0; 
    rparams    = [l; sr; los]; 
    % time series duration
    d1         = '20180101';
    d2         = '20190111';
    dn_all     = [datenum(d1, 'yyyymmdd'):12:datenum(d2, 'yyyymmdd')]';
    % time series connectivity
    gidx_orig  = [[1:length(dn_all)-1]' [2:length(dn_all)]'];
    dc_orig    = dn_all(gidx_orig); 
    du         = sort(unique(dc_orig(:))); 
    nvels      = length(du)-1;
    % velocity and DEM error magnitude
    vel        = 0; 
    dz         = 30; 
    Go         = diff(dn_all(gidx_orig)')'*0.00274; 
    
% for loop paremeters
    % dates of dem error introduction
    ccdates    = [du(1)-12; du+1]; 
    nccdates   = length(ccdates); 
    % number of rand baselines generated per DEM error introduction step
    ntests     = 1; 
    % misc
    mb_all1    = {};
    mv_all1    = {};
    mb_all     = {};
    mv_all     = {};
    mbnl_all   = {};
    mvnl_all   = {};
    mtnl_all   = {};
    options    = optimoptions(@lsqnonlin,'Display','off', 'OptimalityTolerance', 1e-14, 'Algorithm', 'levenberg-marquardt');
    timerVal   = tic;
    %tv2        = tic; 



% loop through date of DEM introduction, remake G every time
for i = 1:nccdates
    % for loop parameters
    c      = ccdates(i); 
    btwint = zeros(length(dc_orig), 1); 
    aftint = ones(length(dc_orig), 1); 
    mb1    = [];
    mv1    = []; 
    mb     = [];
    mv     = [];
    mbnl   = [];
    mvnl   = [];
    mtnl   = [];
    
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
    G         = Go(gidx(:,1));  
    %dzi       = make_variable_dz(dz, dc, 0.00274); 
        
    % define "real" deformation
        def   = G*vel; 
        intsr = def;
        
    
     for k =1:ntests
         % generate random imaging geometry
         bli   = [0; (randn(length(dn_all),1)*70)]; 
         bli   = bli(gidx); 
         bli   = bli(:,2)-bli(:,1); 
         % baseline to phase change of DEM error (in meters)
         blgi  = (2*bli)./(rparams(2).*sind(rparams(3))); 
%          lz    = length(find(aftint == 0)); 
%          blgi  = [ones(lz,1)*10e-10; blgi(lz+1:end)]; 
         Gbl   = [G blgi];
        
        % Add noise to ints (if noise is wanted, uncomment " + n;")
        %noise weight: 0.001~0.5cm noise, 0.01~5cm noise
        nw     = 0.001;    
        n      = (randn(length(dc), 1)).*nw;
        intsiz = intsr;% + n;
        
        % add baseline to ints
        iz     = blgi.*dz;
        intsiz(aftint == 0) = intsiz(aftint == 0) + iz(aftint == 0);  % 'lidar'
        %intsiz(aftint == 1) = intsiz(aftint == 1) + iz(aftint == 1); % 'srtm'
        
        % linear inversion
        mest1         = inv(Gbl'*Gbl)*Gbl'*intsiz; 
        
        % iterative linear inversion
        % for srtm
%         X2    = [];
%         for l = 0:length(blgi)-1
%             Gbl      = [G [zeros(l,1); blgi(l+1:end)]]; 
%             mestl    = inv(Gbl'*Gbl)*Gbl'*intsiz; 
%             rl       = (Gbl*mestl)-intsiz; 
%             X2       = [X2; rl'*rl]; 
%         end
%         [mn, midx]   = min(X2); 
%         Gbl          = [G [zeros(midx-1,1); blgi(midx:end)]]; 
%         mest         = inv(Gbl'*Gbl)*Gbl'*intsiz; 
        
        % for lidar
        X2    = [];
        lbg   = length(blgi); 
        for l = fliplr(0:length(blgi)-1)
            Gbl      = [G [blgi(1:lbg-l); zeros(l,1)]]; 
            mestl    = inv(Gbl'*Gbl)*Gbl'*intsiz; 
            rl       = (Gbl*mestl)-intsiz; 
            X2       = [X2; rl'*rl]; 
        end
        [mn, midx]   = min(X2); 
        Gbl          = [G [blgi(1:midx); zeros(length(G)-midx,1)]]; 
        mest         = inv(Gbl'*Gbl)*Gbl'*intsiz; 
        
        
        % nonlinear inversion
        % for srtm
%           mdc1 = mean(dc')-min(dc(:));
%           mdc  = (mdc1-median(mdc1))./0.1; 
%           gss  = ((c - mean(dc(:))))./0.1;
        % for lidar
          mdc1 = mean(dc')-min(dc(:));
          mdc  = (mdc1-median(mdc1));%./0.1; 
          gss  = ((c - mean(dc(:))));%./0.1;
        fun  = @(r) (G*r(1)) + ((-tanh(mdc-r(2))'+1).*blgi*r(3).*0.5) - intsiz;
        x0   = [0 gss 0];
        x    = lsqnonlin(fun, x0, [-Inf -Inf -Inf], [Inf Inf Inf], options);
        %disp([num2str(x0(2)) ' ' num2str(x(2))]); 
        
        % ignore anomolously large velocity estimates
        if mest(1) < 50       
            mb1      = [mb1; mest1(2)]; 
            mv1      = [mv1; mest1(1)]; 
            mb       = [mb; mest(2)]; 
            mv       = [mv; mest(1)]; 
            mbnl     = [mbnl; x(3)]; 
            mvnl     = [mvnl; x(1)];
            mtnl     = [mtnl; x(2)];
        end
     end
     
    % aggregate all velocity and DEM error estimates
    mb_all1   = [mb_all1; mb1]; 
    mv_all1   = [mv_all1; mv1]; 
    mb_all   = [mb_all; mb]; 
    mv_all   = [mv_all; mv]; 
    mbnl_all = [mbnl_all; mbnl]; 
    mvnl_all = [mvnl_all; mvnl];
    mtnl_all = [mtnl_all; mtnl];
    
    % approx time estimate for script to run
    if find(i == 1:10:nccdates)
        elapsedTime = toc(timerVal);
        percentDone = (i./nccdates);
        timeLeft    = (elapsedTime/percentDone)-elapsedTime; 
        str         = datestr(seconds(timeLeft),'HH:MM:SS');
        disp([num2str(round(percentDone*100)) '% done,  time left: ' str]); 
    end
        
   
end

% el2 = toc(tv2);
% disp([num2str(round((datenum(d2,'yyyymmdd')-datenum(d1,'yyyymmdd'))./365)) ' years, ' num2str(el2) 'sec']);

%% Plot

m0    = [];
s1    = [];
mbl1  = [];
sbl1  = [];
m     = [];
s     = [];
mbl   = [];
sbl   = [];
mnl   = [];
snl   = [];
mblnl = [];
sblnl = [];
m1    = [];
m2    = [];
m3    = [];
b1    = [];
b2    = [];
b3    = [];
for j = 1:nccdates
    mvi1 = cell2mat(mv_all1(j))*100; % *100 for m -> cm
    m0   = [m0; mean(mvi1)']; 
    s1   = [s1; std(mvi1)']; 
    mbl1 = [mbl1; mean(cell2mat(mb_all1(j)))]; 
    sbl1 = [sbl1; std(cell2mat(mb_all1(j)))]; 
    m1   = [m1 mvi1]; 
    b1   = [b1 cell2mat(mb_all1(j))]; 
    
    mvi = cell2mat(mv_all(j))*100; % *100 for m -> cm
    m   = [m; mean(mvi)']; 
    s   = [s; std(mvi)']; 
    mbl = [mbl; mean(cell2mat(mb_all(j)))]; 
    sbl = [sbl; std(cell2mat(mb_all(j)))]; 
    m2   = [m2 mvi]; 
    b2   = [b2 cell2mat(mb_all(j))]; 

    
    mvinl = cell2mat(mvnl_all(j))*100; % *100 for m -> cm
    mnl   = [mnl; mean(mvinl)']; 
    snl   = [snl; std(mvinl)']; 
    mblnl = [mblnl; mean(cell2mat(mbnl_all(j)))]; 
    sblnl = [sblnl; std(cell2mat(mbnl_all(j)))]; 
    m3    = [m3 mvinl]; 
    b3    = [b3 cell2mat(mbnl_all(j))]; 
end

% plot mean and standard deviation of velocity
% fig = figure; 
% left_color  = [0 0 0];
% right_color = [0.1 0.5 0.3];
% set(fig,'defaultAxesColorOrder', [left_color; right_color]);
% hold on; box on;
% yyaxis left
% plot([du(1)-1000 du(end)+1000], [0 0], '--', 'color', [0.6 0.6 0.6]); 
% errorbar(ccdates, m, s, 'linewidth', 1);
% plot(ccdates, m, '.k', 'markersize', 10); 
% ylabel('Mean Velocity & std (cm/yr)'); 
% yyaxis right
% errorbar(ccdates+1, mbl, sbl, 'linewidth', 1, 'color', [0.1 0.5 0.3]);
% plot(ccdates+1, mbl, '.', 'markersize', 10, 'color', [0.1 0.5 0.3]); 
% ylim([-40 40]); 
% xlabel('Date of DEM introduction'); 
% ylabel('Mean DEM error & std (m)'); 
% datetick; 
% xlim([du(1)-100 du(end)+100]); 

% close all; 

% plot nonlin vs lin
yl = 1; 
figure('units', 'normalized', 'outerposition', [.3 .3 .5 .4]); hold on; box on;

subplot(1,2,1); 
% figure; 
hold on; box on;
plot([du(1)-100 du(end)+100], [0 0], '--', 'color', [0.6 0.6 0.6]); 
errorbar(ccdates, m0, s1, 'linewidth', 1, 'color', [0.0 0.0 0.0]);
plot(ccdates, m0, '.k', 'markersize', 10, 'color', [0.0 0.0 0.0]); 
errorbar(ccdates+2, m, s, 'linewidth', 1, 'color', [0.1 0.7 0.3]);
plot(ccdates+2, m, '.k', 'markersize', 10, 'color', [0.1 0.7 0.3]); 
errorbar(ccdates-2, mnl, snl, 'linewidth', 1, 'color', [0.1 0.3 0.8]);
plot(ccdates-2, mnl, '.', 'markersize', 10, 'color', [0.1 0.3 0.8]); 
ylabel('Mean Velocity & std (cm/yr)'); 
xlabel('Date of DEM introduction'); 
datetick('x', 'yyyy', 'keeplimits'); 
xlim([ccdates(1)-100 ccdates(end)+100]); 
ylim([-yl yl])


subplot(1,2,2);
% figure; 
hold on; box on;
plot([du(1)-100 du(end)+100], [0 0], '--', 'color', [0.6 0.6 0.6]); 
errorbar(ccdates, mbl1, sbl1, 'linewidth', 1, 'color', [0.0 0.0 0.0]);
plot(ccdates, mbl1, '.k', 'markersize', 10, 'color', [0.0 0.0 0.0]); 
errorbar(ccdates, mbl, sbl, 'linewidth', 1, 'color', [0.1 0.7 0.3]);
plot(ccdates, mbl, '.k', 'markersize', 10, 'color', [0.1 0.7 0.3]); 
errorbar(ccdates, mblnl, sblnl, 'linewidth', 1, 'color', [0.1 0.3 0.8]);
plot(ccdates, mblnl, '.', 'markersize', 10, 'color', [0.1 0.3 0.8]); 
xlabel('Date of DEM introduction'); 
ylabel('DEM error mean & std (m)'); 
datetick('x', 'yyyy', 'keeplimits'); 
xlim([ccdates(1)-100 ccdates(end)+100]); 

ylim([-5 35]);










