
%keyboard

close all
clear

sr         = 9.48e5; 
los        = 43.0; 
l          = 0.055; 
dt         = 12; 
dy         = 12/365; %decimal years per interferogram
dz         = 30; 
bpr        = 70; %bp range
nints      = 32; 
true       = [0;dz]; %true model
timerVal   = tic;
ntrials    = 1; 
t111 = 0; t222 = 0; t333 = 0; 

t4 = tic; 

for j = 1:ntrials

  Bl0   = randn(nints+1,1)*bpr;
  Bl0   = diff(Bl0);
  Bl    = (4*pi*Bl0)./(l*sr*sind(los));
  G     = [ones(nints, 1)*dy*4*pi/l Bl];  %for vel term-  dimensional analysis:  radians/int = meters/year  * radians/meter  * years/int
  ints  = G*true;
  
    for i = 1:nints 
        use            = [1:i-1 i+1:nints];
        Gi             = G(use,:);
        intsi          = ints;
        intsi(i+1:end) = 0;
        intsi          = intsi(use);
        
        % linear
        t1 = tic; 
        Ginv           = inv(Gi'*Gi)*Gi'; 
        m              = Ginv*intsi; 
        Vavg_all1(i,j) = m(1)*100;  
        Havg_all1(i,j) = m(2);
        t11  = toc(t1); 
        t111 = t111 + t11; 
        
        
        % grid search
        t2 = tic; 
        Gii = Gi; 
        X2 = [];
        for k = -1:length(use)-1
            if k ~= -1
                Gii(end-k:end,2) = 0; 
            end
            mestl    = inv(Gii'*Gii)*Gii'*intsi; 
            rl       = (Gii*mestl)-intsi; 
            X2       = [X2; rl'*rl]; 
        end 
        [mn, midx]   = min(X2); 
        midx         = midx-1; 

        Gii          = [Gi(:,1) [Gi(1:end-midx,2); zeros(midx,1)]]; 
        m2           = inv(Gii'*Gii)*Gii'*intsi; 
        Vavg_all2(i,j) = m2(1)*100; 
        Havg_all2(i,j) = m2(2); 
        t22  = toc(t2); 
        t222 = t222 + t22; 
        
        
        % lev-mar 
        t3 = tic; 
        dd         = [1:length(use)]; 
        mdc        = [(dd-min(dd))];
        mdc        = (mdc - median(mdc))*100; 
        gss        = (-0.5+i-mean(dd))*100;
        fun        = @(r) (Gi(:,1)*r(1)) + ((-tanh(mdc-r(2))'+1).*Gi(:,2)*r(3).*0.5) - intsi;
        x0         = [0 gss 0];
        options    = optimoptions(@lsqnonlin,'Display','off', 'OptimalityTolerance', 1e-9); %, 'Algorithm', 'levenberg-marquardt');
        x          = lsqnonlin(fun, x0, [-Inf -Inf -Inf], [Inf Inf Inf], options);
        Vavg_all3(i,j) = x(1)*100; 
        Havg_all3(i,j) = x(3); 
        t33  = toc(t3); 
        t333 = t333 + t33; 
        
    end
    
    if find(j == 1:ntrials/10:ntrials)
        elapsedTime = toc(timerVal);
        percentDone = (j./ntrials);
        timeLeft    = (elapsedTime/percentDone)-elapsedTime; 
        str         = datestr(seconds(timeLeft),'HH:MM:SS');
        disp([num2str(round(percentDone*100)) '% done,  time left: ' str]); 
    end
end
toc(timerVal)
t444 = toc(t4); 

%% Plot
close all; 
d1         = datenum('20180101', 'yyyymmdd');
d2         = d1+(12*(nints-1)); 
x          = d1:12:d2; 
h          = figure; 
lc         = [0.1 0.1 0.1];
rc         = [0 0.5 0];
set(h,'defaultAxesColorOrder',[lc; rc]);
set(gca, 'FontName', 'Arial', 'fontsize', 12);
hold on; box on; 
yyaxis left
set(gca, 'FontName', 'Arial', 'fontsize', 12);
hold on; 
options = struct('x_axis', x, 'error', 'std', 'handle', h, 'color_area', ...
    [0.5 0.5 0.5], 'color_line', [0.5 0.5 0.5], 'alpha', 0.2, 'line_width', 1, ...
    'markertype', 'x', 'markersize', 7);
    plot_areaerrorbar(Vavg_all1', options); hold on; 
options = struct('x_axis', x, 'error', 'std', 'handle', h, 'color_area', ...
    [0.8 0.8 0.8], 'color_line', [0.4 0.4 0.4], 'alpha', 0.2, 'line_width', 1, ...
    'markertype', '.', 'markersize', 10);
    plot_areaerrorbar(Vavg_all2', options); hold on; 
options = struct('x_axis', x, 'error', 'std', 'handle', h, 'color_area', ...
    [0.3 0.3 0.3], 'color_line', [0.3 0.3 0.3], 'alpha', 0.2, 'line_width', 1, ...
    'markertype', 'o', 'markersize', 8);
    plot_areaerrorbar(Vavg_all3', options); hold on; 
ylim([-0.5 0.5])
ylabel('Inferred Average Velocity (cm/year)'); 
xlabel('Date of DEM Error Removal'); 
yyaxis right
set(gca, 'FontName', 'Arial', 'fontsize', 12);
hold on
options = struct('x_axis', x, 'error', 'std', 'handle', h, 'color_area', ...
    [0 0.5 0], 'color_line', [0 0.6 0], 'alpha', 0.2, 'line_width', 1, ...
    'markertype', 'x', 'markersize', 7);
    plot_areaerrorbar(Havg_all1', options); hold on; 
options = struct('x_axis', x, 'error', 'std', 'handle', h, 'color_area', ...
    [0 0.8 0], 'color_line', [0 0.5 0], 'alpha', 0.2, 'line_width', 1, ...
    'markertype', '.', 'markersize', 10);
    plot_areaerrorbar(Havg_all2', options); hold on; 
options = struct('x_axis', x, 'error', 'std', 'handle', h, 'color_area', ...
    [0 0.3 0], 'color_line', [0 0.4 0], 'alpha', 0.2, 'line_width', 1, ...
    'markertype', 'o', 'markersize', 8);
    plot_areaerrorbar(Havg_all3', options); hold on; 
datetick; 
ylim([-18 48])
ylabel('Inferred DEM Error (m)'); 
xlim([x(1)-24 x(end)+24]); 
yticks([-15 0 15 30 45])
yl = get(gca, 'ylabel');
ylp = get(yl, 'Position');
ext=get(yl,'Extent');
set(yl, 'rotation', 270, 'VerticalAlignment','middle','Position',ylp+[ext(3)-15 0 0]);
xticklabels({'2018', ' ', ' ', ' ', ' ', ' ', ' ', '', '', '', '', '', '2019'})











































