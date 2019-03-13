% p222f870
load(['baselines.mat']); 
load(['meancor_bl_dates_area2_HH.mat']); 
l = 0.236; sr = 8.5e5; los = 38.7; 

dn_all = baselines2.dn_all; 
bl_all = baselines2.bl_all; 
d      = meancor_bl_dates; 
gidx   = d.good_cor_idx; 
dc     = d.dateCombos; 
bl     = abs(d.bl); % have to make bl positive for inversion to work
dta    = diff(dn_all); 

bl     = bl(gidx,:); 
dc     = dc(gidx,:); 
dcv    = dc(:);
du     = sort(unique(dcv));
nints  = length(dc); 


% load data exported from ts_v1.m in packrat ~/matlab/
% int
% vr = 'invproj_v1'; 
% load(['phs_box_vec_' vr '.mat']); 
% phs_alla = phs_box_vec.phs_all; 
% phs_all = phs_alla(106:145, 153:194, :); % needs to be y:y, x:x, :
% single point
load('defv1.mat'); % def
phs_all = reshape(def, 1, 1, 32); 



% Build G matrix
Gz   = zeros(length(dc), 2);
G    = zeros(length(dc), 1);
xa   = [];
ya   = [];
for i = 1:length(dc)
    dd       = dc(i,2) - dc(i,1); 
    G(i)     = dd;  
    Gz(i,1)  = dd; 
    Gz(i,2)  = (4.*pi.*bl(i))./(l.*sr.*sind(los)); 
    xa       = [xa; find(du == dc(i,1))];
    ya       = [ya; find(du == dc(i,2))];
end
d2y     = 0.00274; % convert days to years
Gz(:,1) = Gz(:,1).*d2y; 
G       = G.*d2y; 
Gzi = inv(Gz'*Gz)*Gz'; 
Gi  = inv(G'*G)*G'; 

nvel = max(ya); 

% Build G for nonlin inversion 
% make symbolic variables
blg    = [(4.*pi.*bl)./(l.*sr.*sind(los))]; 
z = sym('z'); 
v = sym('v');
svar  = [v; z];

% create residual function
G_m_ = double(G)*v;

% covariance matrix
Cd   = eye(nints).*(1.^2); 
Cdi  = double(inv(sqrt(Cd))); 


% do inversion pixel-by-pixel
npx = size(phs_all, 2); 
npy = size(phs_all, 1); 
mz = []; 
m  = [];
zz = [];
mznl = [];
lamda = 1; 
W = eye(2).*lamda; 
for q = 1:npy
    for o = 1:npx
    vect    = [phs_all(q,o,:)];
    vect    = vect(:); 
% loop through Gs
r_all = [];
mnlvel = [];
mnlz = [];
for l = 0 %1:nvel
    for j = 1:length(dc)
        if xa(j) > l
            G_m_(j) = G_m_(j)+blg(j)*z; 
        end
    end

    % residual function
    f = Cdi*(G_m_ - vect); 

    % create Jacobian 
    J = [];
    for i = 1:2
        vii = svar(i); 
        Ji  = diff(f, vii);
        J   = [J Ji]; 
    end

    % eps
    var0  = zeros(2, 1); 
    %var0  = 1:nvels+1; 
    ep    = 1e-8; 

    % LM method
    [varf, k, Cm, X2]   = LMLSQ_project(f, var0, J, ep, svar); 
    mnlvel = [mnlvel; varf(1:end-1)]; 
    mnlz   = [mnlz; varf(end)]; 

    r_all = [r_all; X2];
    %disp(l); 
end 

[mnnnn, idx] = min(r_all);
mznl(q,o) = mnlvel(idx); 


% linear inversion 
        mestz   = Gzi*vect;
        mest    = Gi*vect;
        mz(q,o) = mestz(1); 
        m(q,o)  = mest(1); 
        zz(q,o) = mestz(2); 
    end
    disp(num2str(round((npx*(q-1)+o)/(npy*npx), 2)));
end


% figure; hold on; 
% bar([m, mz, mznl]); 

%% plot 
%close all; 
figure('units', 'normalized', 'outerposition', [.9 .9 .3 .9]); hold on;
mnn = -30; 
mxx = 30; 
subplot(3,1,1); hold on; 
pcolor((m)); shading flat; colorbar; hold on; axis equal; 
%plot([zs(3) zs(3) zs(4) zs(4) zs(3)], [zs(1) zs(2) zs(2) zs(1) zs(1)], 'k-'); 
caxis([mnn mxx]); 
title('no dem correction'); 
subplot(3,1,2); 
pcolor((mz)); shading flat; colorbar; hold on; axis equal; 
%plot([zs(3) zs(3) zs(4) zs(4) zs(3)], [zs(1) zs(2) zs(2) zs(1) zs(1)], 'k-'); 
title('dem correction'); 
caxis([mnn mxx]); 
subplot(3,1,3); 
pcolor((mznl)); shading flat; colorbar; hold on; axis equal; 
title('non linear')
caxis([mnn mxx]); 

linkaxes
axis([0 npx 0 npy]);
colormap bluewhitered
























