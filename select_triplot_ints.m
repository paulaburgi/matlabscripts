% select_triplot_ints
clear
%close all

pf_fol  = '/data/pmb229/isce/p222f870/'; 
cintfol = 'mostcombos/';
%cintfol = 'HVcombos/';
datafol = [pf_fol 'data/']; 
load([datafol 'analysis/meancor_bl_dates_HH.mat']); 
    cor = meancor_bl_dates.meancor; 
    bl  = meancor_bl_dates.bl; 
    dcombos = meancor_bl_dates.dateCombos; 

cor_lim = 0.20; 
good_cor_idx = find(cor >= cor_lim); 
good_cor     = cor(good_cor_idx); 
good_bl      = bl(good_cor_idx); 
good_dcombos = dcombos(good_cor_idx,:); 

bad_cor_idx  = find(cor < cor_lim); 
bad_cor     = cor(bad_cor_idx); 
bad_bl      = bl(bad_cor_idx); 
bad_dcombos = dcombos(bad_cor_idx,:); 

figure; 
subplot(1,2,1); 
hist(good_bl); 
xlabel('spatial baselines'); 
subplot(1,2,2); 
hist(diff(good_dcombos')); 
xlabel('temporal baselines'); 