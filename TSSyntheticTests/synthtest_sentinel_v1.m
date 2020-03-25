% TopsDEMTest.m 
% This script creates synthetic InSAR data that contains a DEM error for
% different portions of the full time series. 


%% Variable Definitions
sr         = 9.48e5;   % Slant Range (m)
los        = 43.0;     % Line-of-sight angle (degrees)
l          = 0.055;    % Wavelength (m)
r2m        = (4*pi)/l; % Radians to meters conversion
dt         = 12;       % Temporal spacing (days)
dy         = dt/365;   % Decimal years per interferogram
dz         = 30;       % True DEM error (m) 
bpr        = 70;       % Baseline range (m)
nints      = 30;       % Number of interferograms
nw         = 0;      % Noise weight (cm)
ntrials    = 1000;     % Number of iterations
true       = [0;dz];   % True model


%% Synthetic test
% Iterate over number of trials
Vavg_all   = zeros(nints,ntrials); 
Havg_all   = zeros(nints,ntrials); 
for j = 1:ntrials
    
    % Generate and scale random baseline distribution
    Bl0   = randn(nints+1,1)*bpr;
    Bl0   = diff(Bl0);
    Bl    = (r2m*Bl0)./(sr*sind(los));
    
    % Create design matrix
    G     = [ones(nints, 1)*dy*r2m Bl]; 
    
    % True data + noise
    ints  = G*true;
    noise = diff(randn(nints+1,1)*nw);
    
    % Iterate over number of interferograms, removing DEM error at each
    % time step
    for i = 1:nints 
        % Do not use interferograms that span date of DEM removal
        use            = [1:i-1 i+1:nints];
        Gi             = G(use,:);
        intsi          = ints;
        % Set interferograms after date of DEM error removal to zero
        intsi(i+1:end) = 0;
        intsi          = intsi(use) + noise(use);

        % Invert for average velocity and DEM error
%         Ginv          = inv(Gi'*Gi)*Gi'; 
        Ginv          = (Gi'*Gi)\Gi'; 
        m             = Ginv*intsi; 
        Vavg_all(i,j) = m(1)*100;  
        Havg_all(i,j) = m(2);
    end
end

Vm = mean(Vavg_all,2)';
Ve = std(Vavg_all,0,2)';
Hm = mean(Havg_all,2)';
He = std(Havg_all,0,2)';
% a = [a; nw Vm(1)+Ve(1)]; 
% keyboard; 


%% Plot
% Define date interval
d1         = datenum('20180101', 'yyyymmdd');
d2         = d1+(dt*(nints-1)); 
x          = d1:dt:d2; 
xp         = [x, fliplr(x)]; 

% Create figure, assign axis colors
h          = figure; 
lc         = [0.5 0.5 0.5];
rc         = [0 0.5 0];
set(h,'defaultAxesColorOrder',[lc; rc]);
set(gca, 'FontName', 'Arial', 'fontsize', 12)
hold on; box on; 
% Plot Vavg mean and standard deviation
yyaxis left  
    fill(xp, [Vm+Ve,fliplr(Vm-Ve)], lc, 'edgecolor', lc, 'LineStyle', '-', 'FaceAlpha', 0.33);
    plot(x, Vm, '.-', 'color', lc, 'markersize', 10);
    ylabel('Inferred Average Velocity (cm/year)'); 
    %ylim([-0.53 0.53]); 
% Plot Havg mean and standard deviation
yyaxis right 
    fill(xp, [Hm+He,fliplr(Hm-He)], rc, 'edgecolor', rc, 'LineStyle', '-', 'FaceAlpha', 0.33);
    plot(x, Hm, '.-', 'color', rc, 'markersize', 10);
    ylabel('Inferred DEM Error (m)');
    ylim([-9 38]); 
xlabel('Date of DEM Error Removal (years)'); 
xlim([d1-24 d2+24]);
datetick('x','keeplimits');













% %% Plot with noise
% 
% close all; 
% 
% figure; 
% set(gca, 'FontName', 'Arial', 'fontsize', 12)
% hold on; box on; 
% % Plot Vavg mean and standard deviation
% lc2 = [0.2 0.4 0.6];
% fill([x(1) x(end) x(end) x(1)], [Vm(1)+Ve(1), Vm(end)+Ve(end),Vm(end)-Ve(end), Vm(1)-Ve(1)], lc2, 'edgecolor', lc2, 'LineStyle', '-', 'FaceAlpha', 0.33);    
% fill(xp, [Vm+Ve,fliplr(Vm-Ve)], lc, 'edgecolor', lc, 'LineStyle', '-', 'FaceAlpha', 0.33);
% plot(x, Vm, '.-', 'color', lc, 'markersize', 10);
% ylabel('Inferred Average Velocity (cm/year)'); 
% ylim([-1.8 1.8]); 
% xlabel('Date of DEM Error Removal (years)'); 
% xlim([d1-24 d2+24]);
% datetick('x','keeplimits');





































