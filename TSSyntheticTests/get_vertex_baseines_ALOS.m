

%clear 
ccdir = '/data/pmb229/other/clearcuttingTStest/'; 
cd(ccdir); 


% % Query Vertex API, finding all ALOS SAR frames in Cascadia
% global (CAREFUL, takes forever) ALSO this is only for 2010-2011... 
%     [x,y]=system(['curl https://api.daac.asf.alaska.edu/services/search/param?polygon=' ... 
     %            '-124.4,40.8,-122.0,40.6,-121.4,48.8,-124.8,48.8,-124.4,40.8\&platform=' ...
%                 'ALOS\&processingLevel=L1.0\&offNadirAngle=34.3\&output=csv']);
% cascadia     
     %[x,y]=system(['curl https://api.daac.asf.alaska.edu/services/search/param?polygon=' ... 
     %            '-124.4,40.8,-122.0,40.6,-121.4,48.8,-124.8,48.8,-124.4,40.8\&platform=' ...
     %            'ALOS\&processingLevel=L1.0\&offNadirAngle=34.3\&output=csv']);

     %save('api_result_ALOS_global.mat', 'y'); % 'api_result_ALOS_cascadia.mat','api_result_ALOS_global.mat'

% load query results, if already ran
%     load('api_result_ALOS_global.mat'); 
%     yg = y; 
    load('api_result_ALOS_cascadia.mat'); 
    yc = y; 
    
%     y = yg; 

% 31 variables, 30 commas seperating them. 
c   = strfind(y,','); % find commas
nc  = length(c);      % number of commas
nv  = 30;             % number of variables (30 FOR CASCADIA, 31 FOR GLOBAL) 
m   = 1:nc;           % comma indexes
nf = (nc/nv)-1;       % number of frames

% element of interest (Bp=29, path=6, frame=7, date=8)
vb  = 29;
vp  = 6; 
vf  = 7; 
vd  = 8; 
vl1 = 15; 
vl2 = 16; 
vl3 = 17; 
vl4 = 18; 
vl5 = 19; 
vl6 = 20; 
vl7 = 21; 
vl8 = 22; 

% loop through each from and get info
bl = []; p = []; f = []; d = []; l = [];
for i = 1:nf
    mi    = m(i); 
    % baseline
    t     = y(c(vb+(mi*nv)) : c(vb+1+(mi*nv))); 
    bl    = [bl;  str2double(t(3:end-2))]; 
    % path
    t     = y(c(vp+(mi*nv)) : c(vp+1+(mi*nv))); 
    p     = [p;   str2double(t(3:end-2))]; 
    % frame
    t     = y(c(vf+(mi*nv)) : c(vf+1+(mi*nv))); 
    f     = [f;   str2double(t(3:end-2))]; 
    % date
    t     = y(c(vd+(mi*nv)) : c(vd+1+(mi*nv))); 
    d     = [d;   datenum(t(3:12))]; 
    % lat/lon
    t1    = y(c(vl1+(mi*nv)) : c(vl1+1+(mi*nv))); 
    t2    = y(c(vl2+(mi*nv)) : c(vl2+1+(mi*nv))); 
    t3    = y(c(vl3+(mi*nv)) : c(vl3+1+(mi*nv))); 
    t4    = y(c(vl4+(mi*nv)) : c(vl4+1+(mi*nv))); 
    t5    = y(c(vl5+(mi*nv)) : c(vl5+1+(mi*nv))); 
    t6    = y(c(vl6+(mi*nv)) : c(vl6+1+(mi*nv)));
    t7    = y(c(vl7+(mi*nv)) : c(vl7+1+(mi*nv)));
    t8    = y(c(vl8+(mi*nv)) : c(vl8+1+(mi*nv)));
    l     = [l;  str2double(t1(3:end-2)) str2double(t2(3:end-2)) str2double(t3(3:end-2)) ...
             str2double(t4(3:end-2)) str2double(t7(3:end-2)) str2double(t8(3:end-2)) ...
             str2double(t5(3:end-2)) str2double(t6(3:end-2)) str2double(t1(3:end-2)) str2double(t2(3:end-2))]; 
    
end

% find unique path/frame combos
pf          = [p f];
[C,ia,ic]   = unique(pf, 'rows'); 
% find number of each unique path/frame combo
npf         = accumarray(ic(:),1); 
% only look at path/frames with more than 10 frames
gidx = find(npf > 15); % 15

%% Plot
    figure; hold on; 
    xlim([datenum('2007-01-01', 'yyyy-mm-dd') datenum('2012-01-01', 'yyyy-mm-dd') ]); 
    ylim([-1.5e4 1e4]);
    dcmObj = datacursormode;
set(dcmObj, 'UpdateFcn', @datacursorprecision);

d_bl   = {}; 
pf_all = [];
lall   = [];
for i= 1:length(C)
    if ismember(i,gidx)
        % find indices for each unique p/f combo
            idxi = find(ic == i); 
        % pull out baseline and dates for all indices for each p/f
            [di,didx]  = sort(d(idxi,:)); 
            bli        = bl(idxi,:); 
            bli        = bli(didx,:); 
            pfi        = pf(idxi,:); 
            pfi        = pfi(didx,:); 
            li         = l(idxi,:); 
            li         = li(1,:); 
        % plot them
            plot(di, bli, 'k.', 'markersize', 20); 
            %pause(0.1); 
            oops
            plot(di, bli, 'k.', 'color', [0.8 0.8 0.8], 'markersize', 5); 
        % concatenate all date/baselines for each p/f combo
            d_bl   = [d_bl; di bli];
            pf_all = [pf_all; pfi(1,:)];
            lall   = [lall; li]; 
    end
end


%% Make baseline combos

dt_lim = 150; % days
bl_lim = 1000; % meters

gidx_all = {};
for i = 1:length(d_bl)
    dbli          = cell2mat(d_bl(i)); 
    [dn_all,sidx] = sort(dbli(:,1), 'ascend'); 
    bl_all        = dbli(:,2);
    bl_all        = bl_all(sidx); 
    
    dcom          = combnk(dn_all,2);
    bcom          = combnk(bl_all,2);

    diffd         = abs(diff(dcom'));
    bc            = find(diffd<dt_lim);
    dcom_tb       = dcom(bc,:); 
    bcom_tb       = bcom(bc,:); 

    diffb         = abs(diff(bcom_tb')); 
    dc            = find(diffb<bl_lim);
    gcombos       = dcom_tb(dc,:);
    [~,gidx]      = ismember(gcombos, dn_all); 
    gidx_all      = [gidx_all; gidx];
end

% save(['ALOS_all_params_d_bl_d' num2str(dt_lim) '-bl' num2str(bl_lim) '.mat'], 'd_bl'); 
% save(['ALOS_all_params_gidx_all_d' num2str(dt_lim) '-bl' num2str(bl_lim) '.mat'], 'gidx_all'); 




%% plot for andes seminar

close all; 

figure; hold on; box on; 
plot([datenum('01012006', 'ddmmyyyy') datenum('01052012', 'ddmmyyyy')], [0 0], 'k--'); 


for i =1:length(d_bl)
    dbli = cell2mat(d_bl(i)); 
    gidx = cell2mat(gidx_all(i)); 
    dni  = dbli(:,1); 
    dni  = dni(gidx); 
    bli  = dbli(:,2); 
    bli  = bli(gidx); 
    for j = 1:length(bli)
        plot(dni(j,:), bli(j,:), 'color', [0.7 0.7 0.7]); %, cmap(i,:)); 
    end
    plot(dni, bli, '.', 'markersize', 10, 'color', [0.7 0.7 0.7]); %, cmap(i,:));
end
dbli = cell2mat(d_bl(57)); 
gidx = cell2mat(gidx_all(57)); 
dni  = dbli(:,1); 
dni  = dni(gidx); 
bli  = dbli(:,2); 
bli  = bli(gidx); 
for j = 1:length(bli)
    plot(dni(j,:), bli(j,:), 'k', 'linewidth', 2); 
end
plot(dni, bli, 'k.', 'markersize', 20);
datetick; 
xlabel('Date'); 
ylabel('Perp Baseline (m)'); 
xlim([datenum('01012007', 'ddmmyyyy') datenum('01052011', 'ddmmyyyy')]); 
ylim([-6500 6500]);





%% Lat/lon 

%close all
% states = shaperead('states', 'UseGeoCoords', true);
% sn    = [2 9 12 23 25];
states = shaperead('usastatehi', 'UseGeoCoords', true);
sn     = [5 12 28 37 47];
s = states([sn]);

figure; hold on; box on; 
f=fill([-125.8 -125.8 -120.5 -120.5], [40 49.5 49.5 40], [0.85 0.9 0.95]);
geoshow(s, 'DefaultFaceColor', [0.8 0.8 0.8], 'DefaultEdgeColor', [0.4 0.4 0.4]);
plot(lall(:,2:2:end)', lall(:,1:2:end)', 'color', [0.0 0.0 0.0], 'linewidth', 1)
axis equal
axis([-125.8 -120.5 40 49.5]); 
%set(gca, 'color', [0.85 0.9 0.95]);
title('script in "get vertex baselines ALOS.m"'); 
























    








