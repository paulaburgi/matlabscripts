
% SEE IF I CAN SIMPLIFY THIS SCRIPT TO BARE BONES, SAME RESULT? 

clear

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
    vel        = 0; 
    dz         = 30; 
    Go         = diff(dn_all(gidx_orig)')'*0.00274; 
    
% for loop paremeters
    % dates of dem error introduction
    ccdates    = [du(1)-12; du+1]; 
    nccdates   = length(ccdates); 
    % number of rand baselines generated per DEM error introduction step
    ntests     = 1000; 
    % misc
    mb_all     = {};
    mv_all     = {};
    
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
         bli   = [0; (randn(length(dn_all),1))]; 
         bli   = bli(gidx); 
         bli   = bli(:,2)-bli(:,1); 
         % baseline to phase change of DEM error (in meters)
         blgi  = 70*(2*bli)./(rparams(2).*sind(rparams(3))); 
         Gbl   = [G blgi];
        
        intsiz = intsr;
        
        % add baseline to ints
        iz     = blgi.*dz;
        intsiz(aftint == 0) = intsiz(aftint == 0) + iz(aftint == 0);  % 'lidar'
        
        % linear inversion
        mest     = inv(Gbl'*Gbl)*Gbl'*intsiz; 
        
        
        % ignore anomolously large velocity estimates
        if mest(1) < 50      
            mb       = [mb; mest(2)]; 
            mv       = [mv; mest(1)]; 
        end
     end
     
    % aggregate all velocity and DEM error estimates
    mb_all   = [mb_all; mb]; 
    mv_all   = [mv_all; mv]; 
    
    
end


%% Plot


m     = [];
s     = [];
mbl   = [];
sbl   = [];
m2    = [];
b2    = [];
for j = 1:nccdates
    mvi = cell2mat(mv_all(j))*100; % *100 for m -> cm
    m   = [m; mean(mvi)']; 
    s   = [s; std(mvi)']; 
    mbl = [mbl; mean(cell2mat(mb_all(j)))]; 
    sbl = [sbl; std(cell2mat(mb_all(j)))]; 
    m2  = [m2 mvi]; 
    b2  = [b2 cell2mat(mb_all(j))]; 
end

% plot mean and standard deviation of velocity
fig = figure; 
left_color  = [0 0 0];
right_color = [0.1 0.5 0.3];
set(fig,'defaultAxesColorOrder', [left_color; right_color]);
hold on; box on;
yyaxis left
plot([du(1)-1000 du(end)+1000], [0 0], '--', 'color', [0.6 0.6 0.6]); 
errorbar(ccdates, m, s, 'linewidth', 1);
plot(ccdates, m, '.k', 'markersize', 10); 
ylabel('Mean Velocity & std (cm/yr)'); 
yyaxis right
errorbar(ccdates+1, mbl, sbl, 'linewidth', 1, 'color', [0.1 0.5 0.3]);
plot(ccdates+1, mbl, '.', 'markersize', 10, 'color', [0.1 0.5 0.3]); 
ylim([-40 40]); 
xlabel('Date of DEM introduction'); 
ylabel('Mean DEM error & std (m)'); 
datetick; 
xlim([du(1)-100 du(end)+100]); 



%%

%% test 2

clear
nints = 20; 
G0    = [ones(nints, 1)];
ints  = zeros(nints, 1);

Vavg_all = []; 
for i = 0:nints
    if i == 0
        Gi    = G0; 
    elseif i == nints
        Gi      = G0; 
    elseif i ~= 0 && i ~= nints
        Gi      = G0([1:i-1 i+1:end], :); 
    end
    
    Vavg  = [];
    ints  = zeros(nints, 1);
    for j = 1:1000
        Bl    = randn(nints,2); 
        Bl    = diff(Bl')'; 
        ints(1:i) = Bl(1:i); 
        if i == 0
            intsi = ints; 
            Gi    = [Gi Bl]; 
        elseif i == nints
            intsi   = ints; 
            Gi    = [Gi Bl]; 
        elseif i ~= 0 && i ~= nints
            intsi   = ints([1:i-1 i+1:end]); 
            Gi      = [Gi Bl([1:i-1 i+1:end])]; 
        end
    
        Ginv = inv(Gi'*Gi)*Gi'; 
        m    = Ginv*intsi; 
        if m(1) < 90e7
            Vavg = [Vavg; m(1)]; 
        else
            Vavg = [Vavg; 0]; 
        end
    
    end
    Vavg_all = [Vavg_all Vavg]; 
end


close all; 
figure; hold on; box on; 
plot(0:20, zeros(21,1), 'color', [0.7 0.7 0.7], 'linewidth', 2); 
Vm = mean(Vavg_all,1)';
Vs = std(Vavg_all);
errorbar(Vm, Vs, 'k'); 



