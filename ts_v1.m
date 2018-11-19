% p222f870
load(['/data/pmb229/isce/p222f870/data/baselines/baselines.mat']); 
load(['/data/pmb229/isce/p222f870/data/analysis/meancor_bl_dates_area2_HH.mat']); 
l = 0.236; sr = 8.5e5; los = 38.7; 
pf_fol  = '/data/pmb229/isce/p222f870/mostcombos/';
cd(pf_fol); 

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

% large region (v2)
y1 = 1e3;
y2 = 1.9e3;
x1 = 6.5e2;
x2 = 2.2e3; 
zs = [733 748 841 853]; 
v = 'v2'; 
% small region (v1)
% y1 = 975; 
% y2 = 1224; 
% x1 = 1193; 
% x2 = 1365; 
% v  = 'v1'; 
% region cut in 2006 with little errors
x1 = 1600; 
x2 = 2000; 
y1 = 700;
y2 = 900;
zs = [78 82 148 162];
zs = [105 109 181 186];
% region cut during ts
x1 = 1300; 
x2 = 1600; 
y1 = 700;
y2 = 900; 
zs = [118 121 133 143]; % yyxx
% zs = [108 108 143 143]; % 1 pixel in cleared region
% zs = [127 127 181 181]; % 1 pixel in stable region
v = 'invproj_v1'; 

% size
npy = (y2-y1)+1; 
npx = (x2-x1)+1;  

za = []; 

% if exist([pf_fol 'timeseries/phs_box_vec_' v '.mat'], 'file')
%     load([pf_fol 'timeseries/phs_box_vec_' v '.mat']); 
%     phs_all = phs_box_vec.phs_all; 
%     npix    = phs_box_vec.npix; 
% else
phs_all = zeros(npy, npx, nints);
% load all data
for i =  1:nints %:nints %round(linspace(1,nints, 16)) %1:nints
        d1 = datestr(dc(i,1), 'yymmdd'); 
        d2 = datestr(dc(i,2), 'yymmdd'); 
        intdir = (['int_' d1 '_' d2 '/']);
 
        % get width and length of cor file
       if i == 1
            x=importdata([intdir 'topophase.cor.geo.vrt']);
            l1 = x{1}; 
            qf = strfind(l1, '"'); 
            nx = str2num(l1(qf(1)+1:qf(2)-1)); 
            ny = str2num(l1(qf(3)+1:qf(4)-1)); 
       end
        
        %intfile
        im = sqrt(-1); 
        filename2 = [intdir 'filt_topophase.unw.geo']; 
        h         = fopen(filename2,'r');
        [F,count] = fread(h,2*nx*ny,'float32');
        status    = fclose(h); 
        rmg       = reshape(F,2*nx,ny); 
        phs       = rmg((nx+1):(nx*2),:);
        mag       = flipud(rmg(1:nx,:)');
        phs       = flipud(phs');
        phsbox    = phs(y1:y2, x1:x2); 
        z         = mean(mean(phsbox(zs(1):zs(2), zs(3):zs(4)))); 
        za = [za; z];
        phs_all(:,:,i) = phsbox-z; 
        disp(i); 
        p = phsbox(:); 
        npix        = length(p); 
        
        if i == nints
            phs_box_vec = struct('y1y2x1x2', [y1 y2 x1 x2]', ...
                          'phs_all', phs_all, 'npix', npix); 
            save([pf_fol 'timeseries/phs_box_vec_' v '.mat'], 'phs_box_vec'); 
        end
end
% end


% build G matrix
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

% do inversion pixel-by-pixel
mz = []; 
m  = [];
zz = [];
lamda = 1; 
for i = 1:npy
    for j = 1:npx
        vect    = [phs_all(i,j,:)];
        vect    = vect(:); 
        mestz   = inv((Gz'*Gz)+eye(2).*lamda)*Gz'*vect;  
        mest    = inv((G'*G))*G'*vect;  
        mzvel   = mestz(1); 
        mvel    = mest(1); 
        mz(i,j) = mzvel; 
        m(i,j)  = mvel; 
        zz(i,j) = mestz(2); 
    end
end

% load first topophase.flat.geo 
% filename2   = ['int_070713_071013/filt_topophase.flat.geo']; %OR .FLAT.GEO 
% fid         = fopen(filename2,'r','native');
% [rmg2,cnt]  = fread(fid,[nx*2,ny],'real*4');
% status      = fclose(fid);
% real        = flipud((rmg2(1:2:nx*2,1:ny))');
% imag        = flipud((rmg2(2:2:nx*2,1:ny))');
% phs         = (angle(real+im*imag));
% mag         = abs(real+im*imag);
% phsbox1     = phs(y1:y2, x1:x2); 
% magbox1     = mag(y1:y2, x1:x2); 

%% plot 
%close all; 
figure; hold on; 
subplot(2,1,1); hold on; 
pcolor((m)); shading flat; colorbar; hold on; axis equal; 
plot([zs(3) zs(3) zs(4) zs(4) zs(3)], [zs(1) zs(2) zs(2) zs(1) zs(1)], 'k-'); 
caxis([-10 10]); 
title('no dem correction'); 
subplot(2,1,2); 
pcolor((mz)); shading flat; colorbar; hold on; axis equal; 
plot([zs(3) zs(3) zs(4) zs(4) zs(3)], [zs(1) zs(2) zs(2) zs(1) zs(1)], 'k-'); 
title('dem correction'); 
caxis([-10 10]); 


linkaxes
axis([0 npx 0 npy]);
























