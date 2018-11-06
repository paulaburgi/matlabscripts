% p222f870
load(['/data/pmb229/isce/p222f870/data/baselines/baselines.mat']); 
load(['/data/pmb229/isce/p222f870/data/analysis/meancor_bl_dates_area2_HH.mat']); 
l = 0.236; sr = 8.5e5; los = 38.7; 

dn_all = baselines2.dn_all; 
bl_all = baselines2.bl_all; 
d      = meancor_bl_dates; 
gidx   = d.good_cor_idx; 
dc     = d.dateCombos; 
bl     = abs(d.bl); % have to make bl positive for inversion to work
dta    = diff(dn_all); 

