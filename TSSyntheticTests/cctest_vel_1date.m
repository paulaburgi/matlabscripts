%% Load data, build basic G

% close all
clear 

ccdir = '/data/pmb229/other/clearcuttingTStest/'; 

% load baseline data & int data
    load([ccdir 'ALOS_p222f870_params_tscombos_old.mat']); 
    l = 0.236; sr = 8.5e5; los = 38.7; 

% extract useful info
    dn_all = params.dn_all; 
    bl_all = params.bl_all; 
    gidx   = params.intcombos; 
    dta    = diff(dn_all); 

% select which ints to work with
    bl     = bl_all(gidx); 
    bl     = bl(:,2)-bl(:,1);
    dc     = dn_all(gidx); 
    du     = sort(unique(dc(:))); 
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
    rdi = find(sum(G)' == 0);
    nrdi = length(rdi);
    G  = [Ga; ones(nrdi, nvels)*(1/(nvels-1))];
    rng = 1:nrdi; 
    for i = rng
        idx = i-1; 
        G(end-idx, rdi(i)) = -1; 
    end
    zrdi = zeros(nrdi, 1); % need to add to end of future vectors

%% Compare data to inv with and without DEM Error
% noise weight
nw  = 1; % nw=1 ~ +/-1cm noise (disp = lamda*phase / 4*pi) 
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
iz     = blg.*dz;   % addition to phase due to topo error
intsrz = intsr + iz; 
intsnz = intsrz + n; 

% calc Gsvd with baseline term 
r2      = rank(Gbl); 
[U,S,V] = svd(Gbl); 
Gsvd    = V(:,1:r2)*inv(S(1:r2,1:r2))*U(:,1:r2)'; 
mzsvd   = Gsvd*intsnz; 
mvel0   = mzsvd(1:end-1); 
mz0     = mzsvd(end); 

%% Add it DEM error with each time step 
% date of clear cutting
    % c = du(1)-1      -> no DEM error in TS
    % c = du(1)+1      -> DEM error present for whole TS
    % c = du(ndates)+1 -> no DEM error in TS
c  = du(11); 

% Initiate ints where only some have additional topo error phase
intsiz = intsr;  
btwint = zeros(length(dc), 1); 
% loop through ints, if DEM error is present after date c, then 
for j = 1:length(dc)
    idx1 = dc(j,1) <= c;
    idx2 = dc(j,2) >= c;
    % add phase to ints that include dates greater than j
    if idx1 && idx2
        btwint(j,1) = 1; 
        intsiz(j) = intsiz(j)+iz(j); 
    elseif ~idx1
        intsiz(j) = intsiz(j)+iz(j); 
    end
end

% ALSO: make it so you take out dates between which cutting happens

mest     = Gsvd*intsiz;
mbl      = mest(end); 
mvel     = mest(1:end-1); 
mz       = mest(end); 


%% nonlinear inversion 

% make symbolic variables
z = sym('z'); 
v = [];
for i = 1:nvels
    if i < 10
        vi = ['v0' num2str(i)]; 
    else
        vi = ['v' num2str(i)]; 
    end
    vi = sym(vi); 
    v = [v; vi]; 
end
svar  = [v; z];

% create residual function
G_m_ = double(G)*v;

% covariance matrix
Cd   = eye(nints+nrdi).*(nw^2); 
Cdi  = double(inv(sqrt(Cd))); 

% loop through Gs
mnlvel = [];
mnlz = [];
r_all = [];
for l = 1%1:nvels
    for j = 1:length(dc)
        if xa(j) > c
            G_m_(j) = G_m_(j)+blg(j)*z; 
        end
    end

    % residual function
    f = Cdi*(G_m_ - double(intsiz)); 

    % create Jacobian 
    J = [];
    for i = 1:nvels+1
        if i == nvels+1
            Ji  = diff(f, z); 
        else
            vii = v(i); 
            Ji  = diff(f, vii);
        end
        J   = [J Ji]; 
    end

    % eps
    %var0  = zeros(nvels+1, 1); 
    var0  = 1:nvels+1; 
    ep    = 1e-3; 

    % LM method
    [varf, k, Cm, X2]   = LMLSQ_project(f, var0, J, ep, svar); 
    f2 = matlabFunction(f);
    %[varf2]   = lsqnonlin(f2, var0); 
    mnlvel = [mnlvel varf(1:end-1)]; 
    mnlz   = [mnlz; varf(end)]; 

    r_all = [r_all; X2];
    disp(l); 
end 

disp(c); 
[m, idx] = min(r_all)
mnlveli = mnlvel(:,idx); 
mnlzi = mnlz(idx); 

% figure; hold on; 
% plot(r_all); 
% plot(idx, m, '*');

%% plot

figure; hold on; box on; 
dd     = diff(du); 
vdates = du(1:end-1)+dd;  
plot(vdates, zeros(nvels, 1), 'k', 'HandleVisibility','off'); 
%plot(vdates, [mvel0], 'linewidth', 2);
plot(vdates, [mvel./100], 'linewidth', 2, 'HandleVisibility','off'); 
plot(vdates, [mvel./100], 'linewidth', 2); 
plot(vdates, [mnlveli/100], 'linewidth', 2); 
yl = ylim; 
plot([c; c], [-100 100], 'k--'); 
plot(vdates(rdi), varf(rdi), '*'); 
xlabel('Date'); 
ylabel('Velocity'); 
datetick('x'); 
xlim([vdates(1)-100 vdates(end)+100]);
ylim([-1 1]);
% ylim(yl); 
p = 3; % precision
legend(['dem error (lin), z=' num2str(mz, p) 'm'], ...
       ['dem error (non-lin), z=' num2str(mnlzi, p) 'm'], ... 
       'missing dates', ...
       'location', 'northeast'); 



% ['no error (lin), z=' num2str(mz0, p) 'm'], ...
% 
% ['dem error start, z=' num2str(dz, p) 'm'], ... 






















    
    
    
    
    
    
    
    
