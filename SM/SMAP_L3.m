cd('/data/pmb229/other/SM/smapL3/'); 
clear 

d    = dir('SMAP_L3*.h5'); 
d    = {d.name}'; 
info = h5info(cell2mat(d(1)));

% load data
for i =1:length(d)
    di     = cell2mat(d(i)); 
    hsm    = hdf5read(di, '/Soil_Moisture_Retrieval_Data_AM/soil_moisture');
    hsmp   = hdf5read(di, '/Soil_Moisture_Retrieval_Data_PM/soil_moisture_pm');
    
    sm_am(:,:,i) = hsm(615:646, 110:152); 
    sm_pm(:,:,i) = hsmp(615:646, 110:152); 
    dt(i,:)      = di(14:21);
end
sm_am(find(sm_am == -9999)) = NaN; 
sm_pm(find(sm_pm == -9999)) = NaN; 
dt                          = datenum(dt, 'yyyymmdd');
small = nan(size(sm_am, 1), size(sm_am, 2), size(sm_am,3)*2);
small(:,:,1:2:end) = sm_am; 
small(:,:,2:2:end) = sm_pm; 
small(small == -9999) = nan; 
dtn = nan(size(dt,1)*2,1); 
dtn(1:2:end) = dt; 
dtn(2:2:end) = dt; 

% idx of point of interest
% p           = [55.105, 19.315];
hlonc       = hdf5read(cell2mat(d(1)), '/Soil_Moisture_Retrieval_Data_AM/longitude');
hlatc       = hdf5read(cell2mat(d(1)), '/Soil_Moisture_Retrieval_Data_AM/latitude');
hlonc = hlonc(615:646, 110:152); 
hlatc = hlatc(615:646, 110:152); 
hlonc(hlonc == -9999) = nan; 
hlatc(hlatc == -9999) = nan; 
% ll          = [hlonc(:) hlatc(:)]; 
% [idx, dist] = dsearchn(ll, p); % check dist is <0.35

% %% plot point of interest
% close all; 
% [i,j] = ind2sub(size(hlonc), idx); 
% y_am  = sm_am(:,i,j); 
% y_pm  = sm_pm(:,i,j); 
% 
% figure; hold on; grid minor; box on; 
% plot(dt, y_am, 'r.', 'markersize', 10); 
% plot(dt, y_pm, 'b.', 'markersize', 10); 
% ylim([min([y_am; y_pm])-0.01  max([y_am; y_pm])+0.01]); 
% xlim([min(dt)-10              max(dt)+10]);
% datetick; 
% xlabel('date'); 
% ylabel('soil moisture'); 
% legend('am', 'pm');


%%
% event 1 timing
dtr = datenum('2018/05/26'); % date of rain
d1  = find(dtn == datenum('05/27/2018')); 
d2  = find(dtn == datenum('09/30/2018')); 
t   = [d1(1); d2(end)];            % may 28 - june 30
dt  = dtn(t(1):t(2)); 
dt1 = dt - dtr; 

ft            = fittype( 'a*exp(-b*x) + c', 'independent', 'x', 'dependent', 'y' );
myopt         = fitoptions( 'Method', 'NonlinearLeastSquares','TolFun',1e-2);
myopt.Display = 'Off';
myopt.Lower   = [0 0 0];
myopt.StartPoint = [dt1(1) 0.1 0];

mags            = nan(size(small,1), size(small,2));
times           = mags;
maglow          = mags;
maghig          = mags;
timerr          = mags;
count            = mags;
for i = 1:size(small, 1)
    for j = 1:size(small, 2)
        d     = squeeze(small(i,j,:));
        d     = d(t(1):t(2)); 
        isf   = isfinite(d); 
        d     = d(isf); 
        dti   = dt1(isf); 
        
        if sum(d == -9999) ~= length(d)
            if length(d) > 10
            [fitresult]      = fit(dti, d, ft, myopt );
            results          = coeffvalues(fitresult);
            confs            = confint(fitresult); % fit confidence

%                 if (sum(confs(:)<0)==0)
                    mags(i,j)   = results(1); %exp(-results(1));
                    times(i,j)  = 1./results(2);
                    %maglow(i,j) = exp(-confs(2,1));
                    %maghig(i,j) = exp(-confs(1,1));
                    %timerr(i,j) = -diff(1./confs(:,2));
%                 else
%                     maglow(i,j) = -10; %threw out because errors included neg.
%                 end
            else
                maglow(i,j) = -30; %threw out because too few pts
            end
        else
            maglow(i,j) = -20; %threw out because too few pts
        end
    end
    
end


%% plot 

bbox = [14.4 27 49.7 60.9]; % latmin latmax, lonmin, lonmax
load('/data/pmb229/other/clearcuttingTStest/WorldHiVectors.mat'); 
a            = [0; find(isnan(lat) == 1)];
in           = find(inpolygon(lat, lon, bbox(1:2), bbox(3:4)) == 1);
%[long, latg] = meshgrid(ulon, ulat); 
long = hlonc; 
latg = hlatc; 
sarbox       = [54.744, 20.887; 57.076, 21.171; 57.859, 17.723; 55.493, 17.336; 54.743, 20.887];

%close all;
figure('units', 'normalized', 'outerposition', [.1 .7 .7 .5]);
subplot(1,2,1); hold on; box on; 
surf(long(:,1:end-5), latg(:,1:end-5), mags(:,1:end-5)); shading flat; view(2); 
plot3(lon, lat, ones(length(lat),1)*100, 'color', 0.5*[1 1 1]); 
% for i = 1:length(a)-1
%     fill(lon(a(i)+1:a(i+1)-1), lat(a(i)+1:a(i+1)-1), 'r'); %0.5*[1 1 1]); 
% end
plot([bbox(3) bbox(3) bbox(4) bbox(4) bbox(3)], [bbox(1) bbox(2) bbox(2) bbox(1) bbox(1)], 'k'); 
plot3(sarbox(:,1), sarbox(:,2), ones(length(sarbox),1)*100, 'k'); 
ylim(bbox(1:2));
xlim(bbox(3:4)); 
title('mags'); 
colorbar; 
caxis([0 0.5]); 

subplot(1,2,2); hold on; box on; 
surf(long(:,1:end-5), latg(:,1:end-5), times(:,1:end-5)); shading flat; view(2); 
plot3(lon, lat, ones(length(lat),1)*100, 'color', 0.5*[1 1 1]); 
plot([bbox(3) bbox(3) bbox(4) bbox(4) bbox(3)], [bbox(1) bbox(2) bbox(2) bbox(1) bbox(1)], 'k'); 
plot3(sarbox(:,1), sarbox(:,2), ones(length(sarbox),1)*100, 'k'); 
ylim(bbox(1:2));
xlim(bbox(3:4)); 
title('times'); 
colorbar; 
%caxis([0 10]);
caxis([0 25]);
linkaxes; 
















