% Sentinel_SynthTest.m
% Written by Paula Burgi, Cornell University, Jan. 2020

% This scipt generates synthetic Sentinel-1 data that simulates the phase
% delay due to forest disturbance during the time series. We invert for
% average velocity assuming a constant DEM error term, and plot the average 
% results for 10,000 trials. 


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
nw         = 0;        % Noise weight (cm)
ntrials    = 10000;    % Number of iterations
true       = [0;dz];   % True model: [avg vel; dem error]

% Create noise constraint design matrix
w         = 1e-4; 
du        = nints+1;
G1        = circshift(eye(nints, du),-1) + eye(nints, du)*-1; 
G1(end,1) = 0; G1(end, end) = 1; 
Gs        = [eye(du).*w zeros(du, 2)];
    
%% Synthetic test

% Define empty model parameter matrices
Vavg_all   = zeros(nints,ntrials); 
Havg_all   = zeros(nints,ntrials); 

% iterate over number of trials
for j = 1:ntrials
    
    % Generate and scale random baseline distribution
    B       = randn(nints+1,1)*bpr;
    Ba(:,j) = B; 
    Bl0     = diff(B);
    Bl      = (r2m*Bl0)./(sr*sind(los));
    
    % Create model parameter design matrix
    G2      = [ones(nints, 1)*dy*r2m Bl]; 
    
    % Augment noise constraint matrix
    G       = [G1 G2];
    Ga      = [G; Gs];
    
    % True data + noise
    ints    = G2*true;
    noise   = diff(randn(nints+1,1)*nw);
        
    % Iterate over number of interferograms, removing DEM error at each
    % time step
    for i = 1:nints 
        % Do not use interferograms that span date of DEM removal
        use            = [1:i-1 i+1:nints];
        Gi             = Ga([use  nints+1:length(Ga)],:);
        Giz            = G(use,:); 
        intsi          = ints;
        % Set interferograms after date of DEM error removal to zero
        intsi(i+1:end) = 0;
        intsi          = intsi(use) + noise(use);

        % Invert for average velocity and DEM error
        Ginv          = (Gi'*Gi)\Giz'; 
        m             = Ginv*intsi; 
        Vavg_all(i,j) = m(end-1)*100; % convert m -> cm 
        Havg_all(i,j) = m(end);
    end
end

Vm = mean(Vavg_all,2)';
Ve = std(Vavg_all,0,2)';
Hm = mean(Havg_all,2)';
He = std(Havg_all,0,2)';



%% Plot

close all; 

% Define variables
d1         = datenum('20180101', 'yyyymmdd'); % date 1 (arbitrary)
d2         = d1+(dt*(nints-1));               % date 2
x          = d1:dt:d2;                        % data iterval
xp         = [x, fliplr(x)];                  
ex         = 1;                               % idx of sample distribution
lc         = [0.5 0.5 0.5];                   % left axis color
rc         = [0 0.5 0];                       % right axis color

% Create figure, assign axis colors
h          = figure('units', 'normalized', 'outerposition', [.1 .7 .35 .8]);
set(h,'defaultAxesColorOrder',[lc; rc]);
axes('position',[0.13 0.1 0.75 0.5]);
set(gca, 'FontName', 'Arial', 'fontsize', 12)
hold on; box on; 

% Plot Vavg 
yyaxis left  
    fill(xp, [Vm+Ve,fliplr(Vm-Ve)], lc, 'edgecolor', lc, 'LineStyle', '-', 'FaceAlpha', 0.33);
    plot(x, Vm, '.-', 'color', lc, 'markersize', 10);
    ylabel('Inferred Avg Velocity (cm/yr)'); 
    plot(x, Vavg_all(:,ex), '--'); 

% Plot Havg
yyaxis right 
    fill(xp, [Hm+He,fliplr(Hm-He)], rc, 'edgecolor', rc, 'LineStyle', '-', 'FaceAlpha', 0.33);
    plot(x, Hm, '.-', 'color', rc, 'markersize', 10);
    ylabel('Inferred DEM Error (m)');
    plot(x, Havg_all(:,ex), '--'); 
xlabel('Date of Forest Disturbance'); 
xlim([d1-24 d2+24]);
datetick('x','keeplimits');

% baseline plot
axes('position',[0.13 0.61 0.75 0.25]);
set(gca, 'FontName', 'Arial', 'fontsize', 12)
xx = [x-(dt/2) x(end)+(dt/2)];
cc = [1 1 1]*0.3;
plot([d1-24 d2+24], [0 0], 'color', 0.9*[1 1 1]); hold on; 
plot(xx, Ba(:,ex), 's-', 'color', cc, 'MarkerEdgeColor', cc, 'MarkerFaceColor',cc, 'markersize', 4); 
xlim([d1-24 d2+24]);
ylim([-250 250]); 
datetick('x','keeplimits');
xlabel('Date'); 
ylabel('Baseline (m)'); 
set(gca, 'XAxisLocation', 'top')





























