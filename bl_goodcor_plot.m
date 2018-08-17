close all

% cascadia
pf_fol  = '/data/pmb229/isce/p222f870/'; 
cintfol = 'mostcombos/';  pol = 'HH'; 
datafol = [pf_fol 'data/']; 
intfol  = [pf_fol cintfol];

cd(intfol); 
intdirs = dir; 
intdirs = {intdirs.name}; 
load([datafol 'analysis/meancor_bl_dates_area2_' pol '.mat']); 
dc      = meancor_bl_dates.dateCombos; 

gidx    = meancor_bl_dates.good_cor_idx; 

load('/data/pmb229/isce/p222f870/data/baselines/baselines.mat');
dn = baselines2.dn_all; 
bl = baselines2.bl_all; 

dcg = dc(gidx,:); 

dci=dcg(:); 
dcu=unique(dci); 


figure; hold on; box on; 
set(gcf, 'Position', [300, 300, 900, 400]) 
plot(dn, bl, 'k.', 'markersize', 15); 

for i=1:length(gidx)
    bl1i = find(dn == dcg(i,1)); 
    bl2i = find(dn == dcg(i,2)); 
    bl1  = bl(bl1i); 
    bl2  = bl(bl2i); 
    plot(dcg(i,:), [bl1 bl2], 'k'); 
end

exints = [733236 733604; 733558 733604; 733374 733604; 733742 733972; ...
          733788 733972; 733282 733604; 733328 733604];
gints  = [733788 733972];
for i=1:length(exints)
    bl1i = find(dn == exints(i,1)); 
    bl2i = find(dn == exints(i,2)); 
    bl1  = bl(bl1i); 
    bl2  = bl(bl2i); 
    plot(exints(i,:), [bl1 bl2], 'r'); 
end

xlab = datestr(dn, 'mm-dd-yyyy'); 
set(gca,'xtick',dn); 
set(gca,'xticklabel',xlab);
set(gca,'XTickLabelRotation',45);