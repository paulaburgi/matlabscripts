% Synthetic time series for DEM error introduced mid-time series. Baselines
% generated randomly, to simulate Sentinel or Nisar-like imaging
% geometries. 

%close all
%clear 
ccdir = '/data/pmb229/other/clearcuttingTStest/'; 
cd(ccdir); 
keyboard; 

%% Set variables 

% Static variables
    % Satellite parameters
    l          = 0.055; 
    sr         = 9.48e5; 
    los        = 43.0; 
    rparams    = [l; sr; los]; 
    % velocity and DEM error magnitude
    vel        = 0; 
    dz         = 30; 
    % number of rand baselines generated per DEM error introduction step
    ntests     = 5000; 
    % durations to test
    d1         = datenum('20100101', 'yyyymmdd'); 
    durs       = datenum(2010,1:61,1); % 5 years = 61, 10 years = 121
    durs       = durs(1:2:end); 
    durs       = durs(4:end); 
    ms_all     = [];
    mm_all     = [];
    
for a = 1:length(durs)    
        % time series duration
        d2         = durs(a);
        dn_all     = [d1 : 12 : d2]';
        % time series connectivity
        gidx_orig  = [[1:length(dn_all)-1]' [2:length(dn_all)]'];
        dc_orig    = dn_all(gidx_orig); 
        du         = sort(unique(dc_orig(:))); 
        nvels      = length(du)-1;
        Go         = diff(dn_all(gidx_orig)')'*0.00274; 

        % for loop paremeters
        % dates of dem error introduction
        ccdates    = [du(1)-12; du+1]; 
        nccdates   = length(ccdates); 

        % misc
        mb_all     = {};
        mv_all     = {};
        timerVal   = tic;
        m   = [];
        s   = [];
        mbl = [];
        sbl = [];

    %% Loop 

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
             Gbl   = [G blgi];

            % Add noise to ints (if noise is wanted, uncomment " + n;")
            %noise weight: 0.001~0.5cm noise, 0.01~5cm noise
            nw     = 0.01;    
            %n      = (randn(length(dc), 1)).*nw;
            intsiz = intsr;% + n;

            % add baseline to ints
            iz     = blgi.*dz;
            intsiz(aftint == 1) = intsiz(aftint == 1) + iz(aftint == 1); 

            % inversion
            mest     = inv(Gbl'*Gbl)*Gbl'*intsiz; 

            % ignore anomolously large velocity estimates
            if mest(1) < 50       
                mb       = [mb; mest(2)]; 
                mv       = [mv; mest(1)]; 
            end
         end

        % aggregate all velocity and DEM error estimates
        mb_all = [mb_all; mb]; 
        mv_all = [mv_all; mv]; 
        m      = [m; mean(mv)];
        s      = [s; std(mv)]; 
        mbl    = [mbl; mean(mb)];
        sbl    = [sbl; std(mb)]; 

    end

    ms_all = [ms_all;  mean(s(3:end-2))*100];
    mm_all = [ mm_all; min(s(3:end-2))*100 max(s(3:end-2))*100]; 
    tc = num2str(round(toc(timerVal))); 
    nw = now;
    nw = datetime(nw,'ConvertFrom','datenum'); 
    nw = datestr(nw); 
    disp([num2str(a) '/' num2str(length(durs)) ' (duration: ' tc 'sec), time now:' nw(end-8:end)]); 
end

%% Plot

close all; 
figure; hold on; box on; 
hold on; 
d = diff([ones(length(durs),1)*d1 durs']')'/365; 
% plot(d, ms_all, 'b.-', 'markersize', 10); 
errorbar(d, ms_all, ms_all-mm_all(:,1), mm_all(:,2)-ms_all, '.-', 'markersize', 10); 
xlabel('Time series duration (years)'); 
ylabel('Mean standard Deviation (cm/year)'); 
ylim([0 max(ms_all)+0.5]); 
xlim([0 max(d)+0.4]); 












