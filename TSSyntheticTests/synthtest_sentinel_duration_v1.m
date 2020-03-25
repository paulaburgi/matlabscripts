clear

sr         = 9.48e5; 
los        = 43.0; 
l          = 0.055; 
dt         = 12; 
dy         = 12/365; %decimal years per interferogram
dz         = 30; 
Vavg_all   = []; 
bpr        = 70; %bp range
true       = [0;dz]; %true model
nw         = 0; 

nintsi     = 15:5:115; %15:5:155;
%nintsi     = 16:5:31;
d          = [((1+nintsi)*12)/365]';
smm        = [];
smn        = [];
for k = 1:length(nintsi)
    nints     = nintsi(k); 
    Vavg_all  = [];
    
    for j = 1:500

      Bl0   = randn(nints+1,1)*bpr;
      Bl0   = diff(Bl0);
      Bl    = (4*pi*Bl0)./(l*sr*sind(los));
      G     = [ones(nints, 1)*dy*4*pi/l Bl];  %for vel term-  dimensional analysis:  radians/int = meters/year  * radians/meter  * years/int
      ints  = G*true;
      noise = diff(randn(nints+1,1)*nw);

        for i = 1:nints 
            use            = [1:i-1 i+1:nints];
            Gi             = G(use,:);
            intsi          = ints;
            intsi(i+1:end) = 0;
            intsi          = intsi(use) + noise(use);


            Ginv          = inv(Gi'*Gi)*Gi'; 
            m             = Ginv*intsi; 
            Vavg_all(i,j) = m(1)*100;  
            %Havg_all(i,j) = m(2);
        end
    end
    s     = [std(Vavg_all(2:end-1, :)')']; 
    smm  = [smm; min(s) max(s)]; 
    %smn = [smn; mean(s)];
    disp(num2str(k)); 
end

% noise addition
% % figure; hold on; 
% plot(d,smm(:,1), ':', 'color', 'r'); 
% plot(d,smm(:,2), ':', 'color', 'b'); 
% 
% keyboard; 

%% Plot
% close all
%h = figure('units', 'normalized', 'outerposition', [.1 .3 .4 .66]); %hold on; box on; 

line(d, smm(:,1), 'color', 'r'); 
line(d, smm(:,2), 'color', 'b'); 

xlim([0 5.4]);
ylim([0 0.92]);
xlabel('Duration of Time Series (years)');
ax1 = gca; 
set(gca, 'FontName', 'Arial', 'fontsize', 13);
% plot(d(4), smm(4,1), '.', 'color', 'r', 'markersize', 10); hold on; 
% plot(d(4), smm(4,2), '.', 'color', 'b', 'markersize', 10); hold off; 
ax1.XColor = 'k';
ax1.YColor = 'k';
ax1_pos = ax1.Position; 
ax2 = axes('Position', ax1_pos, 'XAxisLocation','top','YAxisLocation','right','Color','none');
d2 = nintsi; 
% line(d2, smm(:,1), 'color', [0.5 0 0]); 
% line(d2, smm(:,2), 'color', [0 0 0.5]); 

ax1.YTickLabel = ''; 
xlim([0 163])
ylim([0 0.92]);
xlabel('Number of Interferograms');
ylabel('Inferred Average Velocity (cm/year)'); 
set(gca, 'FontName', 'Arial', 'fontsize', 13);
yl = get(gca, 'ylabel');
ylp = get(yl, 'Position');
ext=get(yl,'Extent');
set(yl, 'rotation', 270, 'VerticalAlignment','middle','Position',ylp+[ext(3)-0.15 0 0]);










































