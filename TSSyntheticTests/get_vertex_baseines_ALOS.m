

clear 
ccdir = '/data/pmb229/other/clearcuttingTStest/'; 
cd(ccdir); 


% % Query Vertex API, finding all ALOS SAR frames in Cascadia
%     %[x,y]=system(['curl https://api.daac.asf.alaska.edu/services/search/param?polygon=' ... 
%     %             '-124.4,40.8,-122.0,40.6,-121.4,48.8,-124.8,48.8,-124.4,40.8\&platform=' ...
%     %             'ALOS\&processingLevel=L1.0\&offNadirAngle=34.3\&output=csv']);
%     %save('api_result_ALOS_cascadia.mat', 'y'); 
% load query results, if already ran
    load('api_result_ALOS_cascadia.mat'); 

% 31 variables, 30 commas seperating them. 
c   = strfind(y,','); % find commas
nc  = length(c);      % number of commas
nv  = 30;             % number of variables 
m   = 1:nc;           % comma indexes
nf = (nc/nv)-1;       % number of frames

% element of interest (Bp=29, path=6, frame=7, date=8)
vb = 29;
vp = 6; 
vf = 7; 
vd = 8; 

% loop through each from and get info
bl = []; p = []; f = []; d = [];
for i = 2:nf
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
end

% find unique path/frame combos
pf          = [p f];
[C,ia,ic]   = unique(pf, 'rows'); 
% find number of each unique path/frame combo
npf         = accumarray(ic(:),1); 
% only look at path/frames with more than 10 frames
gidx = find(npf > 15); 

%% Plot
    figure; hold on; 
    xlim([datenum('2007-01-01', 'yyyy-mm-dd') datenum('2012-01-01', 'yyyy-mm-dd') ]); 
    ylim([-1.5e4 1e4]);
    dcmObj = datacursormode;
set(dcmObj, 'UpdateFcn', @datacursorprecision);

d_bl = {}; 
pf_all=[];
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
        % plot them
            plot(di, bli, 'k.', 'markersize', 20); 
            %pause(0.1); 
            oops
            plot(di, bli, 'k.', 'color', [0.8 0.8 0.8], 'markersize', 5); 
        % concatenate all date/baselines for each p/f combo
            d_bl = [d_bl; di bli];
            pf_all = [pf_all; pfi(1,:)];
    end
end


%% Make baseline combos

dt_lim = 100; % days
bl_lim = 800; % meters

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




























