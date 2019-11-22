% alos data
load('ALOS_p222f870_params_tscombos_old.mat');

close all; 
f=figure('units', 'normalized', 'outerposition', [.1 .3 .3 .7]); 
    left_color  = [0 0 0];
    right_color = [0.5 0.5 0.5];
    set(f,'defaultAxesColorOrder', [left_color; right_color]);

yl = [-4.200 4.000];
subplot('position', [0.08 0.57 0.8 0.4]); hold on; box on; 
    plot([datenum('01012007', 'ddmmyyyy') datenum('01012012', 'ddmmyyyy')], [0 0], 'color', [0.8 0.8 0.8], 'linewidth', 2);
    dn_all = params.dn_all;
    bl_all = params.bl_all;
    ic     = params.intcombos;
    bl     = bl_all(ic); 
    dn     = dn_all(ic); 
    for i=1:length(ic)
        plot(dn(i,:), bl(i,:)/1e3, 'k')
    end
    plot(dn_all, bl_all/1e3, '.k', 'markersize', 10);
    datetick; 
    ylabel('Baseline (km)'); 
    ylim(yl); 
    xlim([datenum('01022007', 'ddmmyyyy') datenum('01092011', 'ddmmyyyy')]);
    
subplot('position', [0.08 0.1 0.8 0.4]); box on;  
    hold on; 
    yyaxis left
    plot([datenum('01012015', 'ddmmyyyy') datenum('01012020', 'ddmmyyyy')], [0 0], 'color', [0.8 0.8 0.8], 'linewidth', 2);
    dn_all = params2.dn_all;
    bl_all = params2.bl_all;
    ic     = params2.intcombos;
    bl     = bl_all(ic); 
    dn     = dn_all(ic); 
%     for i=1:length(ic)
%         plot(dn(i,:), bl(i,:)/1e3, 'k-')
%     end
%     plot(dn_all, bl_all/1e3, '.k', 'markersize', 10);
    datetick; 
    xlabel('Date'); ylabel('Baseline (km)'); 
    ylim(yl); 
    xlim([datenum('01022015', 'ddmmyyyy') datenum('01092019', 'ddmmyyyy')]);
    
    yyaxis right
    plot(dn_all, bl_all, '.', 'markersize', 10, 'color', [0.6 0.6 0.6]);
    for i=1:length(ic)
        plot(dn(i,:), bl(i,:), '-', 'color', [0.6 0.6 0.6])
    end
    ylabel('Baseline (m)'); 
    yl2 = [-210 200]; 
    ylim(yl2)

    for i=1:length(ic)
        plot(dn(i,:), bl(i,:)/(1e3*yl(1)/yl2(1)), 'k-')
    end
    plot(dn_all, bl_all/(1e3*yl(1)/yl2(1)), '.k', 'markersize', 10);
    
yl = get(gca, 'ylabel');
ylp = get(yl, 'Position');
ext=get(yl,'Extent');
set(yl, 'rotation', 270, 'VerticalAlignment','middle','Position',ylp+[ext(3)+0.9 0 0]);

    
    
    
    
    
    
    
    
    