% baseline plot with coherence lines 
clear 
%close all

cd('/data/pmb229/isce/p222f870/mostcombos/'); 

intdirs = dir('int_*'); 
intdirs = {intdirs.name}; 

meancor_all=[];
bl_all = [];
ds1 = []; 
ds2 = []; 

% get mean correlation in area of interest 
for i=1:length(intdirs)
    intdir = cell2mat(intdirs(i)); 
    cd(intdir); 
    %cd('HHHH'); 

   if ~exist('/data/pmb229/isce/p222f870/data/analysis/meancor_bl_dates_area2_HH.mat'); 

    % find size of cor file 
    x=importdata('topophase.cor.vrt');
    l1 = x{1}; 
    qf = strfind(l1, '"'); 
    nx = str2num(l1(qf(1)+1:qf(2)-1)); 
    ny = str2num(l1(qf(3)+1:qf(4)-1)); 

    % import cor file
    filename = 'topophase.cor'; 
    fid         = fopen(filename,'r','native');
    [rmg,count] = fread(fid,[nx,ny*2],'real*4'); 
    status      = fclose(fid); 
    mag         = flipud((rmg(1:nx,1:2:ny*2))');
    phs         = flipud((rmg(1:nx,2:2:ny*2))');

    % isolate area of interest for correlation
    x1 = (nx/10)*7; x2 = (nx/10)*8;
    y1 = (ny/10)*8; y2 = (ny/10)*9; 
    cor = phs(y1:y2, x1:x2); 
    meancor = mean(mean(cor)); 
    meancor_all = [meancor_all; meancor]; 
    end

    % fine baselines 
    [~,bltxt]=system('more isce.log | grep perp_baseline_top');
    bl = str2double(bltxt(29:end)); 
    bl_all = [bl_all; bl];
    ds1 = [ds1; intdir(5:10)];
    ds2 = [ds2; intdir(12:end)];
    dn1 = datenum(ds1, 'yymmdd');
    dn2 = datenum(ds2, 'yymmdd');
    dn = [dn1 dn2];

    %cd ../..
    cd ..
    
end

% if ~exist('/data/pmb229/isce/p222f870/data/analysis/meancor_all.mat'); 
%     save('/data/pmb229/isce/p222f870/data/analysis/meancor_all.mat', 'meancor_all'); 
% end

% load('/data/pmb229/isce/p222f870/data/analysis/meancor_all_area2.mat'); 


load('/data/pmb229/isce/p222f870/data/baselines/baselines2.mat'); 
load('/data/pmb229/isce/p222f870/data/baselines/ds_bl-5e+06m_dl-5e+06m.mat'); 
ds1 = ds.ds1; 
ds2 = ds.ds2; 
d1_all = baselines2.d1_all;  
d2_all = baselines2.d2_all; 
dn_all = baselines2.dn_all; 
bl_all = baselines2.bl_all; 


% plot
figure; hold on; 
set(gcf, 'Position', [300, 300, 900, 400]) 
bl_lim = 1000; 
dt_lim = 500; % 2 years: 730 
dcom = combnk(dn_all,2);
bcom = combnk(bl_all,2);

bll = length(bl_all); 
clr = jet(100); 
clr=[clr(1:55,:); clr(65:end,:)];
clrlim = 10*.5; %max(10*meancor_all); 
clrmax = (0.1.*ceil(clrlim)); 
fct = length(clr)./ clrmax; 

clr2 = ones(length(clr),3); %[.8 .0 .2]
clr2(:,1) = clr2(:,1).*.8; clr2(:,2) = clr2(:,2).*0;  clr2(:,3) = clr2(:,3).*.2; 
%clr = clr2; 

for i=fliplr(1:length(dn))
    idx = find(dcom(:,1) == dn(i,1) & dcom(:,2) == dn(i,2)); 
    xp = dcom(idx,:); 
    yp = bcom(idx,:); 
    clri = ceil(fct*meancor_all(i)); 
    if clri > length(clr)
        clri = length(clr); 
    end
    plot(xp, yp, 'linewidth', 2, 'color', clr(clri,:)); 
    txt1 = (xp(1)+xp(2))/2; 
    txt2 = (yp(1)+yp(2))/2;
    txt3 = num2str(meancor_all(i));
    %text(txt1, txt2, txt3(1:4)); 
end
colormap(clr); 
h=colorbar;
caxis([0 clrmax]); 
ylabel(h, 'Coherence'); 

xlab = cell(1, bll); 

for i=[1:bll]
    plot(dn_all(i),bl_all(i), 'k.', 'markersize', 30); 
    xlab(i)={sprintf('%s-%s-%s',d2_all(i,1:2), d2_all(i,3:4), d2_all(i,5:6))}; 
end

yl = ylim; 

set(gca,'xtick',dn_all); 
set(gca,'xticklabel',xlab);
set(gca,'XTickLabelRotation',45);
set(gca,'ytick',[roundn(yl(1),3):1e3:roundn(yl(2),3)]); 
xlim([min(dn_all(:))-25 max(dn_all(:))+25]); 
ylim([min(bl_all(:))-400 max(bl_all(:))+400]);
%xtick([1:bll])
%xticklabels(xlab)
%xtickangle(75); 
box on; 
xlabel('Date (yy-mm-dd)'); 
ylabel('Baseline (m)'); 
set(gca, 'fontsize', 14); 
ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';
legd = plot(1,1,'.w'); 
legend([legd],sprintf('Temporal Baseline: %g days \nSpatial Baseline:     %g m ', dt_lim, bl_lim), 'location', 'southeast');ds = struct('ds1', ds1, 'ds2', ds2);
ds_name=sprintf('ds_bl-%gm_dl-%gm.mat', bl_lim, dt_lim); 
ds_exist=exist(ds_name); 

title('ALOS data near Eugene, Oregon'); 


% print(gcf, '/data/pmb229/isce/p222f870/data/analysis/cor_bl_plot.pdf', '-dpdf', '-painters'); 

%% plot area of 
clear 
close all

% intdir = 'int_081130_090115/HHHH/'; 
% tit = '11/30/2008 - 01/15/2009';
intdir = 'int_071013_080113/HHHH/'; 
tit = '10/13/2007 - 01/13/2008'; 

x=importdata([intdir 'topophase.cor.vrt']);
    l1 = x{1}; 
    qf = strfind(l1, '"'); 
    nx = str2num(l1(qf(1)+1:qf(2)-1)); 
    ny = str2num(l1(qf(3)+1:qf(4)-1)); 

% import cor file
    filename = [intdir 'topophase.cor']; 
    fid         = fopen(filename,'r','native');
    [rmg,count] = fread(fid,[nx,ny*2],'real*4'); 
    status      = fclose(fid); 
    magcor         = flipud((rmg(1:nx,1:2:ny*2))');
    phscor         = flipud((rmg(1:nx,2:2:ny*2))');
    
% import interferogram
    im = sqrt(-1); 
    filename2 = [intdir 'topophase.flat']; 
    fid          = fopen(filename2,'r','native');
    [rmg2,count] = fread(fid,[nx*2,ny],'real*4');
    status      = fclose(fid);
    real        = flipud((rmg2(1:2:nx*2,1:ny))');
    imag        = flipud((rmg2(2:2:nx*2,1:ny))');
    mag         = abs(real+im*imag);
    phs         = (angle(real+im*imag));

    f = 10;
    x1 = round(((nx/10)*7)/f); x2 = round(((nx/10)*8)/f);
    y1 = round(((ny/10)*8)/f); y2 = round(((ny/10)*9)/f);
    
% zoom out figure
 figure; hold on; 
    hh=gcf; 
    set(hh, 'Position', [50 50 800 400]); 
    ax1 = subplot(1,2,1); hold on; 
    pcolor(phs(1:f:end,1:f:end)); shading flat; 
    plot([x1 x1 x2 x2 x1], [y1 y2 y2 y1 y1], 'k-', 'linewidth', 1.5); 
    axis([0 size(phscor, 2)/f 0 size(phscor, 1)/f]); 
    box on
    caxis([-3.14 3.14]); 
    jetmod = jet(100); 
    jetmod = jetmod(20:80,:); 
    colormap(ax1, jetmod); 
    title(['Interferogram ']); 
    ylabel('Azimuth'); 
    xlabel('Range'); 
    h=colorbar;
    hL = ylabel(h, 'Phase', 'fontsize', 12);
    %set(hL,'Rotation',90);
    ax = gca;
    ax.YAxis.Exponent = 3;
    ax.XAxis.Exponent = 3;
    yticks([0 1 2 3 4 5].*1e3);

    
    %figure; hold on; 
    hh12 = gcf; 
    ax2 = subplot(1,2,2); hold on; 
    pcolor(phscor(1:f:end,1:f:end)); shading flat; 
    plot([x1 x1 x2 x2 x1], [y1 y2 y2 y1 y1], 'k-', 'linewidth', 1.5); 
    axis([0 size(phscor, 2)/f 0 size(phscor, 1)/f]); 
    colormap(ax2, gray);
    box on
    colorbar
    title(['Coherence ']);
    h=colorbar;
    hL = ylabel(h, 'Coherence');
    ylabel(h, 'Coherence', 'fontsize', 12);
    %set(hL,'Rotation',90);
%     ylabel('azimuth'); 
%     xlabel('range'); 
    ax = gca;
    ax.YAxis.Exponent = 3;
    ax.XAxis.Exponent = 3;
    yticks([0 1 2 3 4 5].*1e3);

% zoom in figure 
  figure; hold on; 
    hh2=gcf; 
     set(hh2, 'Position', [10 10 800 400]); 
    ax1 = subplot(1,2,1); hold on; 
    pcolor(phs(1:f:end,1:f:end)); shading flat; 
    %plot([x1 x1 x2 x2 x1], [y1 y2 y2 y1 y1], 'k-', 'linewidth', 1.5);
    axis([x1 x2 y1 y2]); 
    box on
    caxis([-3.14 3.14]); 
    colormap(ax1, jetmod); 
    title(['Interferogram ']); 
    ylabel('Azimuth'); 
    xlabel('Range'); 
    h=colorbar;
    hL = ylabel(h, 'Phase', 'fontsize', 12);
    %set(hL,'Rotation',90);
    ax = gca;
     ax.YAxis.Exponent = 3;
     ax.XAxis.Exponent = 3;
    yticks([4.3 4.4 4.5 4.6 4.7].*1e3);
    
    %figure; hold on; 
    hh22=gcf; 
    ax2 = subplot(1,2,2); hold on; 
    pcolor(phscor(1:f:end,1:f:end)); shading flat; 
    %plot([x1 x1 x2 x2 x1], [y1 y2 y2 y1 y1], 'k-', 'linewidth', 1.5); 
    axis([x1 x2 y1 y2]); 
    colormap(ax2, gray); 
    box on
    colorbar
    title(['Coherence ']);
    h=colorbar;
    hL = ylabel(h, 'Coherence', 'fontsize', 12);
    %set(hL,'Rotation',90);
    %ylabel('azimuth'); 
    %xlabel('range'); 
    ax = gca;
    ax.YAxis.Exponent = 3;
    ax.XAxis.Exponent = 3;
    yticks([4.3 4.4 4.5 4.6 4.7].*1e3); 
    
    set(hh, 'PaperPositionMode', 'auto'); 
    set(hh, 'PaperOrientation', 'landscape');
    set(hh2, 'PaperPositionMode', 'auto'); 
    set(hh2, 'PaperOrientation', 'landscape'); 
    %print(hh, '/data/pmb229/isce/p222f870/data/analysis/bl_cor_loc2_unfilt.jpg', '-djpeg'); 
    %print(hh2, '/data/pmb229/isce/p222f870/data/analysis/bl_cor_loc2zoom_unfilt.jpg', '-djpeg'); 

% print(hh11, '/data/pmb229/isce/p222f870/data/analysis/bl_cor_loc11.jpg', '-djpeg'); 
% print(hh12, '/data/pmb229/isce/p222f870/data/analysis/bl_cor_loc12.jpg', '-djpeg'); 
% print(hh21, '/data/pmb229/isce/p222f870/data/analysis/bl_cor_loc21.jpg', '-djpeg'); 
% print(hh22, '/data/pmb229/isce/p222f870/data/analysis/bl_cor_loc22.jpg', '-djpeg'); 



















    