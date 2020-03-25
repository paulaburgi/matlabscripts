%% Load L2 data all 

cd('/data/pmb229/other/SM/smapL2/'); 
clear

% data directory
d    = dir('SMAP_*.h5'); 
%d    = dir('SMAP_L2*_2018*.h5'); 
d    = {d.name}'; 
info = h5info(cell2mat(d(1)));

% bounding box
bbox = [14.4 27 49.7 60.9]; % latmin latmax, lonmin, lonmax

% name of saved file
%nm   = ['sm_all1.mat'];   % option 1: H-polarization used to find sm
nm   = ['sm_all2.mat'];     % option 2: V-polarization used to find sm
%nm   = ['sm_all3.mat'];   % option 3: Dual-pol

% either load smap data, or load the .mat file containing smap data
latlonv = [];
if ~exist(nm, 'file')
    for i =1:length(d)
        % load each variable
        di     = cell2mat(d(i)); 
        hsm    = hdf5read(di, '/Soil_Moisture_Retrieval_Data/soil_moisture'); % option1 = H-pol (best), 2=v-pol, 3=dual pole
        lati   = hdf5read(di, '/Soil_Moisture_Retrieval_Data/latitude');
        loni   = hdf5read(di, '/Soil_Moisture_Retrieval_Data/longitude');
        cidx   = hdf5read(di, '/Soil_Moisture_Retrieval_Data/EASE_column_index');
        ridx   = hdf5read(di, '/Soil_Moisture_Retrieval_Data/EASE_row_index');
        
        % using row and col indices to reshape data
        tt   = NaN(max(ridx)+1, max(cidx)+1);
        tlat = NaN(max(ridx)+1, max(cidx)+1);
        tlon = NaN(max(ridx)+1, max(cidx)+1);
        for j = 1:length(hsm)
            tt(double(ridx(j))+1, double(cidx(j))+1)   = double(hsm(j));
            tlat(double(ridx(j))+1, double(cidx(j))+1) = double(lati(j));
            tlon(double(ridx(j))+1, double(cidx(j))+1) = double(loni(j));
        end
        
        % determine which data is in the bounding box
        in       = inpolygon(tlat, tlon, bbox(1:2), bbox(3:4)); 
        idx      = find(in == 1); 
        [ii,jj]  = ind2sub(size(tt), idx);
        inhsm    = tt(min(ii):max(ii), min(jj):max(jj)); 
        inlat    = tlat(min(ii):max(ii), min(jj):max(jj)); 
        inlon    = tlon(min(ii):max(ii), min(jj):max(jj)); 
        
        % concatenate everything into a few matrices/vectors
        if ~isempty(tlat(in)) 
            if exist('i2', 'var'); 
                i2 = i2+1; 
            else
                i2 = 1; 
            end
            llin{i2} = cat(3, inhsm, inlat, inlon); 
            dt(i2,:) = di(22:29); 
            inlatv   = inlat(:); 
            inlonv   = inlon(:); 
            latlonv  = [latlonv; inlatv(~isnan(inlatv)) inlonv(~isnan(inlonv))]; 
        end
   
    end
    
    % find unique values of ease grid
    ulat         = unique(latlonv(:,1)); 
    ulon         = unique(latlonv(:,2)); 
    [long, latg] = meshgrid(ulon, ulat); 
    latg         = flipud(latg); 
    hsms         = nan(size(latg)); 
    % place each matrix on same grid
    for i = 1:length(llin)
        lli         = cell2mat(llin(i)); 
        smi            = squeeze(lli(:,:,1)); 
        lti            = squeeze(lli(:,:,2)); 
        lni            = squeeze(lli(:,:,3));  

        [C1, ia1, ib1] = intersect(lti, latg);
        [C2, ia2, ib2] = intersect(lni, long);
        [ii1, ~]     = ind2sub(size(latg), find(latg == C1(end)));
        [~, jj2]     = ind2sub(size(long), find(long == C2(1)));
        
        hsmi           = hsms;
        lati           = hsms; 
        loni           = hsms; 
        
        hsmi(ii1:size(lti, 1)+ii1-1, jj2:size(lti, 2)+jj2-1) = flipud(smi); 
        lati(ii1:size(lti, 1)+ii1-1, jj2:size(lti, 2)+jj2-1) = lti; 
        loni(ii1:size(lti, 1)+ii1-1, jj2:size(lti, 2)+jj2-1) = lni; 
        hsmi(hsmi == -9999) = nan; 
        lati(lati == -9999) = nan; 
        loni(loni == -9999) = nan; 
        small(:,:,i)   = hsmi; 
        latall(:,:,i)  = lati; 
        lonall(:,:,i)  = loni; 
        
    end
    dtn   = datenum(dt, 'yyyymmdd');
    
    
    % Detrend! 
        x          = [1:52]'; 
        d          = day(datetime(dtn, 'convertfrom', 'datenum'),'dayofyear');
        t          = week(datetime(dtn, 'convertfrom', 'datenum'));
        t(t == 53) = 52; 

        % get trimmed means
        for i = 1:52
            idx                  = find(t == i); 
            smalli               = small(:,:,idx); 
            ssp                  = median(smalli, 3, 'omitnan') + 0.02; 
            smalli(smalli > ssp) = NaN; 
            meansm2(:,:,i)       = mean(smalli, 3, 'omitnan');
        end
        m = mean(meansm2, 3, 'omitnan');

        % fit mean
        sm_detrended = small;  
        for i = 1:size(small,1)
            for j = 1:size(small,2)
                y1 = squeeze(meansm2(i,j,:)); 
                if ~isnan(y1)
                    f  = fit(x, y1, 'smoothingspline', 'smoothingparam',  0.01); 
                    a  = f(x); 
                    for k = 1:52
                        idx = find(t == k); 
                        sm_detrended(i,j,idx) = sm_detrended(i,j,idx)-a(k)+m(i,j);
                    end
                end
            end
        end

        %close all; 
        figure; hold on; box on; 
        ii = 16; jj = 13; 
        plot(dtn, squeeze(small(ii,jj,:)), '.'); 
        plot(dtn, squeeze(sm_detrended(ii,jj,:)), '.'); 
        datetick; 
        ylabel('soil moisture'); 
        xlabel('year'); 
        xlim([datenum('2015-01-02') datenum('2020-05-01')]);
        legend('uncorrected', 'corrected', 'location', 'northwest');
        close; 
    
    
    % save as structure
    small = sm_detrended; 
    s     = struct('small', small, 'latall', latall, 'lonall', lonall, ...
                   'ulat', ulat, 'ulon', ulon, 'dtn', dtn); 
    save(nm, 's'); 
else
    load(nm); 
    small  = s.small; 
    smlat  = s.latall; 
    smlon  = s.lonall; 
    dtn    = s.dtn; 
    ulat   = s.ulat; 
    ulon   = s.ulon; 
end




%% inv all

%close all 

% Rain event 1 timing
dtr = datenum('2018/05/26'); % date of rain
d1  = find(dtn == datenum('05/27/2018')); 
d2  = find(dtn == datenum('10/10/2018')); 
t   = [d1(1); d2(end)];            
dt  = dtn(t(1):t(2)); 
dt1 = dt - dtr; 

% Fit options
ft            = fittype( 'a*exp(-b*x) + c', 'independent', 'x', 'dependent', 'y' );
myopt         = fitoptions( 'Method', 'NonlinearLeastSquares','TolFun',1e-2);
myopt.Display = 'Off';
myopt.Lower   = [0 0 0];
myopt.StartPoint = [dt1(1) 0.1 0];

% create variable matrices
mags            = nan(size(small,1), size(small,2));
times           = mags;
maglow          = mags;
maghig          = mags;
timerr          = mags;
count           = mags;
for i = 1:size(small, 1)
    for j = 1:size(small, 2)
        % extract data for each point for the rain date, remove nans. 
        d     = squeeze(small(i,j,:));
        d     = d(t(1):t(2)); 
        isf   = isfinite(d); 
        d     = d(isf); 
        dti   = dt1(isf); 
        
        % fit exponential to each point
        if sum(d == -9999) ~= length(d)
            if length(d) > 10
            [fitresult]      = fit(dti, d, ft, myopt );
            results          = coeffvalues(fitresult);
            confs            = confint(fitresult); % fit confidence

%                 if (sum(confs(:)<0)==0)
                    mags(i,j)   = results(1); 
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

% load international borders
load('/data/pmb229/other/clearcuttingTStest/WorldHiVectors.mat'); 
[long, latg]  = meshgrid(ulon, ulat); 
sarbox       = [54.744, 20.887; 57.076, 21.171; 57.859, 17.723; ...
                 55.493, 17.336; 54.743, 20.887];

%close all;
figure('units', 'normalized', 'outerposition', [.1 .7 .7 .5]);
subplot(1,2,1); hold on; box on; 
surf(long(:,1:end-5), latg(:,1:end-5), mags(:,1:end-5)); shading flat; view(2); 
plot3(lon, lat, ones(length(lat),1)*100, 'color', 0.5*[1 1 1]); 

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


% figure; hold on; 
% small(small == -9999) = nan; 
% surf(long, latg, nansum(small(:,:,274:332), 3)); shading flat; view(2); 
% plot3(lon, lat, ones(length(lat),1)*100, 'color', 0.5*[1 1 1]); 
% plot([bbox(3) bbox(3) bbox(4) bbox(4) bbox(3)], [bbox(1) bbox(2) bbox(2) bbox(1) bbox(1)], 'k'); 
% plot3(sarbox(:,1), sarbox(:,2), ones(length(sarbox),1)*100, 'k'); 
% ylim(bbox(1:2));
% xlim(bbox(3:4)); 




%% Save as mat/geotiff

% .mat
s = struct('long', long(:,1:end), 'latg', latg(:,1:end), ...
    'ulon', ulon(1:end), 'ulat', ulat, ...
    'times', times(:,1:end), 'mags', mags(:,1:end)); 
% option1, 2, or 3
nmat = ['mags_times' nm(end-4) '.mat'];
%save(nmat, 's'); 


% geotiff
% d99 = mags; 
% d99(isnan(d99)) = -99;
% lonlim = [min(long(:)) max(long(:))]; 
% latlim = [min(latg(:)) max(latg(:))]; 
% s = size(d99);
% R = georefcells(latlim, lonlim, s);
% geotiffwrite('mags_2018_all2.tif', d99, R);
% 
% d99 = times; 
% d99(isnan(d99)) = -99;
% lonlim = [min(long(:)) max(long(:))]; 
% latlim = [min(latg(:)) max(latg(:))]; 
% s = size(d99);
% R = georefcells(latlim, lonlim, s);
% geotiffwrite('times_2018_all2.tif', d99, R);











%% show 2 points



% dtr = datenum('2018/05/24'); % date of rain
% d1  = find(dtn == datenum('05/27/2018')); 
% d2  = find(dtn == datenum('10/10/2018')); 
% t   = [d1(1); d2(end)];           
% dt  = dtn(t(1):t(2)); 
% dt1 = dt - dtr; 

% ft            = fittype( 'a*exp(-b*x) + c', 'independent', 'x', 'dependent', 'y' );
% myopt         = fitoptions( 'Method', 'NonlinearLeastSquares','TolFun',1e-2);
% myopt.Display = 'Off';
% myopt.Lower   = [0 0 0];
% myopt.StartPoint = [dt1(1) 0.1 0];

i = 16; j = 13; 
d1             = squeeze(small(i,j,:));
d1             = d1(t(1):t(2)); 
isf            = isfinite(d1); 
d1             = d1(isf); 
dti1           = dt1(isf); 
[fitresult1]   = fit(dti1, d1, ft, myopt );
results1       = coeffvalues(fitresult1);
confs1         = confint(fitresult1); % fit confidence
mags1          = results1(1); %exp(-results(1));
times1         = 1./results1(2);

i = 20; j = 11; 
d2             = squeeze(small(i,j,:));
d2             = d2(t(1):t(2)); 
isf           = isfinite(d2); 
d2             = d2(isf); 
dti2           = dt1(isf); 
[fitresult2]   = fit(dti2, d2, ft, myopt );
results2       = coeffvalues(fitresult2);
confs2         = confint(fitresult2); % fit confidence
mags2          = results2(1); %exp(-results(1));
times2         = 1./results2(2);

figure; hold on; 
plot(dti1, d1, 'b*');
plot(fitresult1, 'b');
plot(dti2, d2, 'r*');
plot(fitresult2, 'r');
plot([times1 times1], [0 0.5], '--b'); 
plot([times2 times2], [0 0.5], '--r'); 

%% Point
% % %% Load L2 data point
% % 
% % cd('/data/pmb229/other/SM/smapL2/'); 
% % clear
% % 
% % d    = dir('SMAP_L2*.h5'); 
% % d    = {d.name}'; 
% % info = h5info(cell2mat(d(1)));
% % 
% % % idx of point of interest
% % %p   = [55.105, 19.315];
% % %p    = [54.4 17.9]; 
% % p    = [55.35, 18.59];
% % 
% % % name of saved file
% % nm   = ['sm_point_' num2str(p(1)) '_' num2str(p(2)) '_option1xxx.mat'];
% % 
% % if ~exist(nm, 'file')
% %     % load data
% %     for i =1:length(d)
% %         di     = cell2mat(d(i)); 
% %         hsm    = hdf5read(di, '/Soil_Moisture_Retrieval_Data/soil_moisture_option1'); % option1 = H-pol (best), 2=v-pol, 3=dual pole
% %         eri    = hdf5read(di, '/Soil_Moisture_Retrieval_Data/soil_moisture_error');         
% %         lati   = hdf5read(di, '/Soil_Moisture_Retrieval_Data/latitude');
% %         loni   = hdf5read(di, '/Soil_Moisture_Retrieval_Data/longitude');
% %          dt(i,:)   = di(22:29);
% % 
% %         ll = [loni lati]; 
% %         [idx, dist] = dsearchn(ll,p); 
% %         if dist < 0.35
% %             [ii,jj] = ind2sub(size(lati), idx); 
% %             y(i)  = hsm(ii,jj); 
% %             er(i) = eri(ii,jj); 
% %         else
% %             y(i)  = NaN; 
% %             er(i) = NaN; 
% %         end
% %     end
% %     dt                = datenum(dt, 'yyyymmdd');
% %     nn = ~isnan(y); 
% %     
% %     % same data point
% %     s = struct('y', y, 'dt', dt, 'nn', nn); % 'er', er); 
% %     save(nm, 's'); 
% % else
% %     load(nm); 
% % end
% % plt = 0; 
% % if plt == 1
% %     fig = figure('units', 'normalized', 'outerposition', [.1 .7 .7 .5]);
% %     hold on;  box on; 
% %     plot(dt(nn), y(nn), 'k.-', 'markersize', 15); 
% %     %errorbar(dt(nn), y(nn), er(nn), 'k.-', 'markersize', 15); 
% %     ylabel('Soil Moisture'); 
% %     xlabel('Date'); 
% %     %xlim([datenum('2018-05-01') datenum('2018-06-15')]);
% %     datetick('x', 'keeplimits');
% %     dcmObj = datacursormode;
% %     set(dcmObj, 'UpdateFcn', @datacursorprecision);
% % end
% % 
% % %% inv point
% % 
% % %close all 
% % 
% % y   = s.y; 
% % dt  = s.dt; 
% % nn  = s.nn; 
% % yn  = y(nn); 
% % dtn = dt(nn); 
% % 
% % % event 1
% % dtr = datenum('2018/05/27'); % date of rain
% % t   = [111; 135]; 
% % y1  = yn(t(1):t(2))'; 
% % dt1 = dtn(t(1):t(2)); 
% % dt1 = dt1 - dtr; 
% % 
% % %y1 = y1-mean(y1(10:end)); 
% % figure; hold on; box on; 
% % plot(dt1, y1, 'o-'); 
% % %datetick('x', 'keeplimits');
% % 
% % 
% % ft            = fittype( 'a*exp(-b*x) + c', 'independent', 'x', 'dependent', 'y' );
% % myopt         = fitoptions( 'Method', 'NonlinearLeastSquares','TolFun',1e-2);
% % myopt.Display = 'Off';
% % myopt.Lower   = [0 0 0];
% % myopt.StartPoint = [dt1(1) 0.1 0];
% % 
% % [fitresult]      = fit(dt1 ,y1, ft, myopt );
% % results          = coeffvalues(fitresult);
% % times            = 1./results(2)
% % plot(fitresult, 'g'); 





























