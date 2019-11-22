close all
clear

sr         = 9.48e5; 
los        = 43.0; 
l          = 0.055; 
dt         = 12; 
dy         = 12/365; %decimal years per interferogram
dz         = 30; 
Vavg_all   = []; 
bpr   = 70; %bp range
nints = 30; %moved this up here, keep things out of loop if possible
true = [0;dz]; %true model

for j = 1:2000

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
        

        Ginv          = inv(Gi'*Gi)*Gi'; 
        m             = Ginv*intsi; 
        Vavg_all(i,j) = m(1)*100;  
        Havg_all(i,j) = m(2);
    end
end

h = figure; hold on; box on; 
options = struct('x_axis', [1:nints]', 'error', 'std', 'handle', gcf, 'color_area', ...
    [0.5 0.5 0.5], 'color_line', [0.5 0.5 0.5], 'alpha', 0.33, 'line_width', 1, ...
    'markertype', '.', 'markersize', 10);
yyaxis left
plot_areaerrorbar(Vavg_all', options);
yyaxis right
hold on
options = struct('x_axis', [1:nints], 'error', 'std', 'handle', gcf, 'color_area', ...
    [0 0.5 0], 'color_line', [0 0.5 0], 'alpha', 0.33, 'line_width', 1, ...
    'markertype', '.', 'markersize', 10);
plot_areaerrorbar(Havg_all', options);
ylim([-40 40])


figure; hold on; box on; 
plot(0:nints+5, zeros(nints+6,1), 'color', [0.7 0.7 0.7], 'linewidth', 2); 
Vm = mean(Vavg_all,2)*100;
Vs = std(Vavg_all')'*100;
errorbar(Vm, Vs, 'k'); 
ylabel('velocity (cm/year)'); 
xlabel('interval of DEM error removal'); 


figure; hold on; box on; 
plot(0:nints+5, dz*ones(nints+6,1), 'color', [0.7 0.7 0.7], 'linewidth', 2); 
Hm = mean(Havg_all,2);
Hs = std(Havg_all')';
errorbar(Hm, Hs, 'k'); 
ylabel('Height'); 
xlabel('interval of DEM error removal'); 












































