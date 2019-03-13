% bl_plot_from_intfols

close all

load('/data/pmb229/isce/p222f870/data/baselines/baselines.mat');
d = dir('int_*'); 
d = {d.name};
ddn = baselines2.dn_all; 
dbl = baselines2.bl_all; 

dn = [];
bl = [];
for i = 1:length(d)
    di = cell2mat(d(i));

    % get dates  a
    d1 = di(5:10); 
    d2 = di(12:end); 
    dn = [dn; datenum(d1, 'yymmdd') datenum(d2, 'yymmdd')]; 
    
    % get baseline for each date
    i1 = round(dbl(find(dn(i,1) == ddn))); 
    i2 = round(dbl(find(dn(i,2) == ddn)));
    bl = [bl; i1 i2]; 
    
end


% dates
da = unique(dn(:));
dc = find(diff(dn')' < 97);
dn2 = dn(dc,:); 
bl2 = bl(dc,:); 

% plot
figure; hold on; 
xlim([min(da)-20 max(da)+20]); ylim([-4.2e3 4.2e3]); 
plot(dn, bl, '.k', 'markersize', 10); 

for i = 1:length(dc)
    % plot
    plot(dn2(i,:), bl2(i,:), 'k'); 
end
datetick('x', 'mm/yyyy');


dc = find(diff(dn')' < 47);
dn2 = dn(dc,:); 
bl2 = bl(dc,:); 
for i = 1:length(dc)
    % plot
    plot(dn2(i,:), bl2(i,:), 'k'); 
end

ik = [da(1:end-1) da(2:end)]; 
ik = [ik(1:7,:); ik(10:12,:); ik(14:end,:); datenum('100605', 'yymmdd') datenum('101206', 'yymmdd')];


close; 
figure('units', 'normalized', 'outerposition', [0 1 1 .5]); hold on;
xlim([min(da)-20 max(da)+20]); ylim([-4.2e3 4.2e3]); 
plot(dn, bl, '.k', 'markersize', 10); 

bldiff = [];
for i=1:length(ik)
    [id, id2] = dsearchn(dn, ik(i,:));
    if id2 == 0
        plot(dn(id,:), bl(id,:), 'b'); 
        bldiff = [bldiff; abs(diff(bl(id,:)))]; 
    else
        %keyboard
    end
end

plot([min(dn(:,1))-100 max(dn(:,1))+100], [0 0], 'k'); 

datetick('x', 'mm/yyyy');


% save variable
% useints= ik;
useints = struct; 
useints.dn = ik; 
useints.bl = bldiff; 
%save('log/useints.mat', 'useints');



































% % get baseline
%     [~, v] = system(['more ' di '/isce.log | grep perp_baseline_top']);
%     bl = [bl; round(str2double(v(strfind(v, '=')+1:end)))];