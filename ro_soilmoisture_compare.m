% on mongoose: /data/rlohman/Sentinel/Saudi/28_67_stack
load cr_forpaula.mat
load output_for_paula.mat

% on mongoose: /data/pmb229/other/temp/sm_compare
dt   = importdata('dt_fromGEE.csv');
sm   = importdata('sm_fromGEE.csv');
sm   = sm.data(:,2);

idx  = find(~isnan(sm)); 
sm   = sm(idx); 
idx2 = find(~isoutlier(sm)); 
sm   = sm(idx2); 
sm   = sm(1:end-6); 

dt2 = [];
for i=2:length(dt)
    t   = cell2mat(dt(i)); 
    sf  = strfind(t, ',');
    d2  = t(sf(1)+1:end-1);
    dt2 = [dt2; d2];
end
dt = datenum(dt2, 'yyyy-mm-dd'); 
dt = dt(idx); 
dt = dt(idx2); 
dt = dt(1:end-6); 
[dt,ud] = unique(dt); 
sm = sm(ud); 

sm = abs(1-sm); 
% m = [1 0; 0 5]; 
% sm2=[];
% for i = 1:length(sm)
%     i2 = m*[1; sm(i)];
%     sm2(i) = i2(2); 
% end



close all; 
fig = figure('units', 'normalized', 'outerposition', [.1 .9 .4 .4]); 
c = 256; 
left_color = [139   0   0]/c;
right_color = [0 128 128]/c;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
hold on; box on; 

yyaxis right
plot(dt, sm, '.-', 'linewidth', 1, 'markersize', 10);
ylabel('Inv. Perpendicular Drought Index (PDI)'); 


yyaxis left
plot(output.dn, cr, '.-', 'linewidth', 1, 'markersize', 10);
xlabel('date'); 
ylabel('InSAR-derived soil moisture');

% plot(output.dn, cr, '.-');
% plot(dt, sm2, '.-');
datetick('x', 'mm-yyyy'); 
xtickangle(45);

orient(fig,'landscape'); 


















