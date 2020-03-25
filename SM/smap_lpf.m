% ***NOTE: a shortened/edited version of this has been transferred to the
%          script: smapinv.m


nm     = 'sm_all2.mat'; 
load(nm); 
small  = s.small; 
dtn    = s.dtn; 

d = day(datetime(dtn, 'convertfrom', 'datenum'),'dayofyear');
t = week(datetime(dtn, 'convertfrom', 'datenum'));
t(t == 53) = 52; 

close all; 
figure; hold on; 
%ii = 16; jj = 13; 
ii = 20; jj = 11; 
% get mins/means/maxes
for i = 1:52; 
    idx   = find(t == i); 
    smalli          = small(:,:,idx); 
    mediansm(:,:,i) = median(smalli, 3, 'omitnan');
    meansm(:,:,i)   = mean(smalli, 3, 'omitnan');
    tmeansm(:,:,i)  = trimmean(smalli, 95,3);
    minsm(:,:,i)    = min(smalli, [],3);

    sd   = std(smalli, 0,3, 'omitnan'); 
    ssp = mediansm(:,:,i) + 0.02; 
    smalli(smalli > ssp) = NaN; 
    meansm2(:,:,i)   = mean(smalli, 3, 'omitnan');
    
    plot(d(idx)./7, squeeze(small(ii,jj,idx)), 'k.'); 
end

% fit mins/means/maxes
sm_detrended = small; 
x1 = [1:52]'; 
for i = 1:size(small,1)
    for j = 1:size(small,2)
        y1 = squeeze(minsm(i,j,:)); 
        %y1 = squeeze(meansm2(i,j,:)); 
        if ~isnan(y1)
            f  = fit(x1, y1, 'smoothingspline', 'smoothingparam',  0.01); 
            a  = f(x1); 
            for k = 1:52
                idx = find(t == k); 
                sm_detrended(i,j,idx) = sm_detrended(i,j,idx)-a(k)-0.03;
            end
        end
    end
end

% plot
plot(1:52, squeeze(minsm(ii,jj,:)), '.-', 'markersize', 15); 
plot(1:52, squeeze(meansm(ii,jj,:)), '.-', 'markersize', 15); 
plot(1:52, squeeze(meansm2(ii,jj,:)), '.-', 'markersize', 15); 
%plot(1:52, squeeze(tmeansm(ii,jj,:)), '.-', 'markersize', 15); 
plot(1:52, squeeze(mediansm(ii,jj,:)), '.-', 'markersize', 15); 
%plot(f); 
box on; 
xlim([0 52.5]); 
%ylim([0.04 0.11]);
ylabel('soil moisture'); 
xlabel('week of the year'); 


% a = f(x1); 
% close all; 
% figure; hold on; 
% for i = 1:52
%     idx   = find(t == i); 
%     small2(:,:,idx) = small(:,:,idx)-a(i);% + 0.08;
%     
% end

%%
%close all; 
figure; hold on; 
ii = 16; jj = 13; 
plot(dtn, squeeze(small(ii,jj,:)), '.'); 
plot(dtn, squeeze(sm_detrended(ii,jj,:)), '.'); 
datetick; 
ylabel('soil moisture'); 
xlabel('year'); 
box on; 
xlim([datenum('2015-01-02') datenum('2020-05-01')]);
legend('uncorrected', 'corrected', 'location', 'northwest');


nm2     = 'sm_all_detrended.mat'; 
s.small = sm_detrended; 
%save(nm2, 's'); 


