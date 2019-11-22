%SM_all_gee2mat

cd /data/pmb229/other/SM/

ps = 'p3'; 

d = importdata(['smap_smos_coh_ndvi_pdi_imerg_' ps '.csv']);

dt = d.textdata(2:end); 
dd = d.data; 
if ps == 'p4'
    dd = [zeros(size(dd,1),1) dd]; 
end
dd(:,1) = dd(:,1)-7; % imerg
dd(:,2) = (dd(:,2)*70)-28; % pdi
dd(:,3) = dd(:,3)+12; %coh
dd(:,4) = (dd(:,4)*15)+21; % ndwi
dd(:,5:6) = dd(:,5:6); % smap/smos

dt2 = datenum(dt); 
df  = find(dt2 == 737182); 
d1  = df(1); 
df  = find(dt2 == 737394); 
d2  = df(end); 
dt2 = dt2(d1:d2,:); 

dd = dd(d1:d2, :); 
dm = repmat(min(dd), size(dd,1),1);
dd = dd-dm; 
dm = repmat(max(dd), size(dd,1),1);
dd = dd./dm; 
t  = repmat([0 0 3 0 1.5 1.5], size(dd,1), 1);
dd = dd+t; 



cm = [0.2 0.2 0.7; ... % imerg
      0.3 0.6 0.3; ...% pdi
      0.4 0.4 0.4; ...% coh
      0.8 0.8 0.1; ...% ndwi snl
      0.5 0.2 0.6; ...% smap
      0.7 0.2 0.2; ...% smos
      ]; 
  
close all; 
f = figure('units', 'normalized', 'outerposition', [.1 .7 .6 .4]); 
hold on; box on;
%dda = {}; 
for i =2:size(dd,2);
    ddi    = dd(:,i);
    n      = ~isnan(ddi); 
    %dda(i) = [dt2(n) ddi(n)]; 
    plot(dt2(n), ddi(n), '.-', 'color', cm(i,:), 'linewidth', 1, 'markersize', 10); 
end

%cygnss data
if ps == 'p1'
    p = [55.105, 19.315]; 
elseif ps == 'p2'
    p = [55.4787, 18.4411]; 
elseif ps == 'p3'
    p = [56.5574, 19.1823];
elseif ps == 'p4'
    p = [56.2374, 19.0494];
end

load('cygnss/cygnssdata.mat'); 
va = []; 
idxa = [];
for i = 121:334
    ln = data(i).lon; 
    lt = data(i).lat; 
    da = data(i).data; 
    da = da(:); 
    [idx,val] = dsearchn([ln(:) lt(:)], p); 
    va = [va; da(idx)]; 
    idxa = [idxa; idx];
end
dts = 737061:737061+364;
dtsi = dts(121:334);

ix = ~isnan(va); 
va = va(ix); 
dtsi = dtsi(ix); 

dm = min(va); 
va = va-dm; 
dm = max(va); 
va = va./dm; 
t  = -1.5; 
va2 = va+t; 

plot(dtsi, va2, '.-', 'color', [0.6 0.3 0.1], 'linewidth', 1, 'markersize', 10); 


set(gca,'YTickLabel',[]);
%ylim([-10 20]); 
datetick; 
xlabel('Date (2018)'); 
ylabel('Soil Moisture'); 
grid minor
ylim([-2 4.5]);

orient(f,'landscape');
%print(gcf, ['plot_' ps '.pdf'], '-dpdf', '-painters');
