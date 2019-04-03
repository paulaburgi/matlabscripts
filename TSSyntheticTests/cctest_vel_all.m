%% Load data, build basic G

 %close all
% clear 

ccdir = '/data/pmb229/other/clearcuttingTStest/'; 

% load baseline data & int data
    load([ccdir 'ALOS_p222f870_params_tscombos.mat']); %17 ints
    % load([ccdir 'ALOS_p222f870_params_tscombos_old.mat']); %32 ints
    %load([ccdir 'ALOS_p222f870_params_tscombos_all.mat']); %210 ints
    l = 0.236; sr = 8.5e5; los = 38.7; 

% extract useful info
    dn_all = params.dn_all; 
    bl_all = params.bl_all; 
    gidx   = params.intcombos; 
    %gidx   = [gidx(1:4,:); gidx(6:end,:)]; 
    %gidx   = [gidx(2:end,:)]; 

% select which ints to work with
    bl     = bl_all(gidx); 
    bl     = bl(:,2)-bl(:,1);
    %bl     = (rand(length(bl), 1)-0.5)*1000; 
    dc     = dn_all(gidx); 
    du     = sort(unique(dc(:))); 
    dta    = diff(du); 
    nvels  = length(du)-1;
    nints  = length(dc); 

% build G matrix
    G  = zeros(nints, nvels);
    xa = [];
    ya = [];
    for i = 1:size(G, 1)
        xi     = dc(i,1); 
        yi     = dc(i,2); 
        x = find(du == xi); 
        y = find(du == yi); 
        
        dtj   = dta(x:y-1);
        for j = 1:length(dtj)
            G(i,x+j-1) = dtj(j); 
        end

        xa = [xa; x];
        ya = [ya; y];
    end
    d2y = 0.00274; % convert days to years
    Ga = G.*d2y; 
    
% Add rows s.t. mean vel = avg vel for each unconstrained date
    rdi    = find(sum(G)' == 0);
    nrdi   = length(rdi);
    avgvel = (1/(nvels-1)); 
    G      = [Ga; ones(nrdi, nvels)*avgvel]; 
    rng    = 1:nrdi; 
    for i = rng
        idx = i-1; 
        G(end-idx, rdi(i)) = -1; 
    end
    zrdi = zeros(nrdi, 1); % need to add to end of future vectors


%% Compare data to inv with and without DEM Error
% noise weight
nw  = 0; % nw=1 ~ +/-1cm noise (disp = lamda*phase / 4*pi) 
n   = [(rand(length(dc), 1)-0.5).*nw; zrdi];

% define time steps and "real" deformation
vel  = ones(nvels, 1).*0; 
def  = Ga*vel; 

% create ints
intsr = [def; zrdi];
ints  = intsr+n;

% add baseline term to G matrix
blg    = [(4.*pi.*bl)./(l.*sr.*sind(los)); zrdi]; 
Gbl    = [G blg]; 
dz     = 30;        % topo error
%dz     = [make_variable_dz(dz, dc, 0.00274); zrdi]; % topo error with tree growth
iz     = blg.*dz;   % addition to phase due to topo error
intsrz = intsr + iz; 
intsnz = intsrz + n; 

% calc Gsvd with baseline term 
r2        = rank(Gbl);
[U,S,V] = svd(Gbl); 
Gsvd    = V(:,1:r2) * inv(S(1:r2,1:r2)) * U(:,1:r2)'; 
mzsvd   = Gsvd*intsrz; 
%test = mzsvd(2:end)-mzsvds; 
mvel0   = mzsvd(1:end-1); 
mz0     = mzsvd(end); 

%% Add it DEM error with each time step 

iz       = blg.*dz; 
mb       = [];
mv       = [];
mb_all   = {};
mv_all   = {};
ccdates  = [du(1)-1; du+1]; 
%ccdates  = [du(1)+1]; 
nccdates = length(ccdates); 
ntests   = 1; 

% loop through date of clear cutting
for i = 1:nccdates
    c = ccdates(i); 
    intsiz = intsr;
    mb = [];
    mv = [];
    btwint = zeros(length(dc), 1); 
    % loops through data values to add Bp affect
    for j = 1:length(dc)
        idx1 = dc(j,1) <= c;
        idx2 = dc(j,2) >= c;
        if idx1 && idx2
            btwint(j) = 1; 
            intsiz(j,1) = intsr(j)+iz(j); 
        elseif idx1
            intsiz(j,1) = intsr(j)+iz(j); 
        else
            intsiz(j,1) = intsr(j); 
        end
    end
    % take out ints that have disturbance between 2 SAR scenes
%     ridx          = find(btwint); 
%     Gblt          = Gbl; 
%     intsiz(ridx)  = 0;
%     Gblt(ridx,:)  = 0; 
%     Gblt(:,ridx)  = 0; 
%     r2            = rank(Gblt);
%     newavgvel = (1/(r2-1)); 
%     Gblt(Gblt == avgvel) = newavgvel; 
%     [U,S,V] = svd(Gblt); 
%     Gsvdt    = V(:,1:r2) * inv(S(1:r2,1:r2)) * U(:,1:r2)'; 
    
    for k = 1:ntests
        n        = [(rand(length(dc), 1)-0.5).*nw; zrdi];
        intsizn  = intsiz + n; 
        %intsizn(ridx) = 0; 
        mest     = Gsvd*intsizn;
        mest     = Gbl\intsizn; 
        mb       = [mb; mest(end)]; 
        mv       = [mv mest(1:end-1)]; 
    end
    mb_all = [mb_all; mb]; 
    mv_all = [mv_all; mv]; 
end


%% Plot

figure; hold on; box on; 
ylim([-20 20]); 
cmaps = jet(nccdates); 
m = [];
for j = 1:nccdates
    mvi = cell2mat(mv_all(j)); 
    for i = 1:ntests
        plot(du, [0; mvi(:,i)], 'color', cmaps(j,:)); 
    end
    m = [m; mean(mean(mvi))];
end
































