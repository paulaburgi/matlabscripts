% pick out phase info from ints

clear
% close all

%% load data
% oregon
     pf_fol  = '/data/pmb229/isce/p222f870/';
     pol = 'HH'; cd([pf_fol 'mostcombos/']);
     %pol = 'HV'; cd([pf_fol 'HVcombos/']);
     datafol = [pf_fol 'data/']; 
    
% sumatra 
%    pf_fol = '/data/pmb229/isce/p446f7190_sumatra/'; 
%    cd([pf_fol 'ints_SRTM']); pol = 'HH'; 
%    datafol = [pf_fol 'data/']; 
    
    load([datafol 'analysis/meancor_bl_dates_area2_' pol '.mat']); 
        datesall  = meancor_bl_dates.dateCombos;     
        idx       = meancor_bl_dates.good_cor_idx;   
        dates     = datesall(idx,:); 
          % if you only want dates on diagonal of triplot, aka daisy chain
               %idx   = meancor_bl_dates.diagonal_idx; 
               %dates = datesall(idx,:); 
        bl        = meancor_bl_dates.bl(idx); 
        datediff  = abs(diff(dates'))'; 
        %diagidx   = meancor_bl_dates.diagonal_idx; 
        %datesdiag = datesall(diagidx, :); 
        
    nints = length(idx); 

    
    
%% Define Areas
    % clear cut
        a  = 1; 
        x1 = 290; x2 = 310;
        y1 = 280; y2 = 290;
%         a  = 2; 
%         x1 = 235; x2 = 255; 
%         y1 = 205; y2 = 215; 
    % forested
        xx1 = 290; xx2 = 310;
        yy1 = 300; yy2 = 310;
%         xx1 = 270; xx2 = 290; 
%         yy1 = 220; yy2 = 230; 

% sumatra 
    % clear cut
       % a=1; 
       % x1 = 320; x2 = 330; 
       % y1 = 125; y2 = 135;  
%         a=2; 
%         x1 = 1144; x2 = 1154; 
%         y1 = 576; y2 = 581; 
    % forested
        %xx1 = 230; xx2 = 240; 
        %yy1 = 105; yy2 = 115;
%         xx1 = 1174; xx2 = 1184; 
%         yy1 = 576; yy2 = 581; 

        
%% For loop
    % for loop variables    
    plotnum=1; 
    pdiff=[];
    error1=[]; 
    cerror=[];
    ferror=[];
    bls=[];
    datediff=[];
    pp1=[];
    pp2=[];
    pe1=[];
    pe2=[]; 
    pdiff2=[];
    pdiff3=[];
    pdiff4=[]; 
    mdiff=[];
    
    blidx=find(bl > 200 & bl < 600)'; 
    figure; 
%         figure; hold on; 
    for i =  1:nints %:nints %round(linspace(1,nints, 16)) %1:nints
        d1 = datestr(dates(i,1), 'yymmdd'); 
        d2 = datestr(dates(i,2), 'yymmdd'); 
        intdir = (['int_' d1 '_' d2 '/']);
 
        % get width and length of cor file
       %if i == 1
            x=importdata([intdir 'topophase.cor.geo.vrt']);
            l1 = x{1}; 
            qf = strfind(l1, '"'); 
            nx = str2num(l1(qf(1)+1:qf(2)-1)); 
            ny = str2num(l1(qf(3)+1:qf(4)-1)); 
        %end
        
        %intfile
         im = sqrt(-1); 
        filename2 = [intdir 'filt_topophase.flat.geo']; %OR .FLAT.GEO 
        fid          = fopen(filename2,'r','native');
        [rmg2,count] = fread(fid,[nx*2,ny],'real*4');
        status      = fclose(fid);
        real        = flipud((rmg2(1:2:nx*2,1:ny))');
        imag        = flipud((rmg2(2:2:nx*2,1:ny))');
        mag         = abs(real+im*imag);
        phs         = (angle(real+im*imag));
        phsbox1         = phs(y1:y2, x1:x2); 
        phsbox2         = phs(yy1:yy2, xx1:xx2); 
        p1 = phsbox1(:); 
        p2 = phsbox2(:); 
            % % figure
            %subplot(3,3,plotnum); 
            plotnum=plotnum+1; 
            %close
%              figure; hold on; %figure('units', 'normalized', 'outerposition', [0 0 1 1]); 
%              pcolor(phs); shading flat; hold on; linkaxes;
%              plot([x1 x1 x2 x2 x1], [y1 y2 y2 y1 y1], 'k-', 'linewidth', 2); 
%              plot([xx1 xx1 xx2 xx2 xx1], [yy1 yy2 yy2 yy1 yy1], 'k-', 'linewidth', 2); 
%             title([d1 '\_' d2 '    B_p=' num2str(round(bl(i))) '    area ' num2str(a)]); 
              %h=colorbar;
%               ylabel(h, 'phase'); 
            %print(gcf, [datafol 'analysis/test' num2str(plotnum-1) '.jpg'], '-djpeg');
            
        % Circular mean
        sinp1 = sum(sin(p1))./length(p1); 
        cosp1 = sum(cos(p1))./length(p1); 
        r1    = sqrt((sinp1.^2) + (cosp1.^2)); 
        sin1  = sinp1/r1; 
        cos1  = cosp1/r1; 
        r11   = sqrt(sin1.^2 + cos1.^2);
        

        
        p1m1  = atan2(sum(sin(p1)), sum(cos(p1))); 

        sinp2 = sum(sin(p2))./length(p2); 
        cosp2 = sum(cos(p2))./length(p2); 
        r2    = sqrt((sinp2.^2) + (cosp2.^2)); 
        sin2  = sinp2/r2; 
        cos2  = cosp2/r2; 
        r22   = sqrt(sin2.^2 + cos2.^2);
        p2m1  = atan2(sum(sin(p2)), sum(cos(p2))); 
        
        p1m   = circ_mean(p1); 
        p2m   = circ_mean(p2); 
        [pd, pd0]    = (circ_dist(p1m, p2m)); 
        pdiff = [pdiff; pd]; 
        re1=r1; 
        re2=r2; 
                
        den = (cos(p2m).^2) + (sin(p2m).^2);
        rl = ((cos(p1m).*cos(p2m))+(sin(p1m).*sin(p2m)))/den;
        img = ((-cos(p1m).*sin(p2m))+(sin(p1m).*cos(p2m)))/den;
        pd2 = atan2(img, rl);
        pdiff2 = [pdiff2; pd2]; 
        
            pp1 = [pp1; p1m1]; 
            pp2 = [pp2; p2m1];
            pe1 = [pe1; p1m];
            pe2 = [pe2; p2m];
        
        % error 
        %e1 = abs(diff([p1m p1e])); 
        %e2 = abs(diff([p2m p2e])); 
        %e1 = 2*pi*(1-r); 
        %e2 = 2*pi*(1-r2); 
        e1  = sqrt(-2.*log(r1));
        e2  = sqrt(-2.*log(r2)); 
        error1 = [error1; sqrt((e1.^2) + (e2.^2))]; 
        cerror = [cerror; e1];
        ferror = [ferror; e2];
        
        % other
        bls = [bls; bl(i)]; 
        datediff = [datediff; abs(dates(i,1)-dates(i,2))];        
      
    end
        
    idx     = find(error1 < 10); 
    bls     = bls(idx, :); 
    pdiff   = pdiff(idx, :); 
    error1  = error1(idx, :); 
    cerror  = cerror(idx, :); 
    ferror  = ferror(idx, :); 
    
    
    
    
%% Plot
    
% plot rms error
    figure; hold on; 
    errorbar(bls, pdiff, error1, '.k', 'markersize', 20); 
    
    ylim([-4 4]); 
    xlabel('Perp Baseline (m)'); 
    ylabel('Phase Difference (radians)'); 
    title(['\Theta_{diff} Between Forested and Cleared Areas (' pol ')']); 
    box on
    xl = xlim; 
    f1 = gcf; 
        
    
    
% plot error individually
    figure; hold on; 
    plot(bls, pdiff, '.k', 'markersize', 20); 
    xlim([xl(1) xl(2)]); 
    c  = round(diff(xl)./350);
    lw = 2; 
    clrs = lines;
    clr1 = clrs(5,:); 
    clr2 = clrs(3,:); 
    oops
    
    for i = 1:length(bls)
        b1 = [bls(i)-c bls(i)-c]; 
        b2 = [bls(i)+c bls(i)+c]; 
        h1 = plot(b2, [pdiff(i)+cerror(i) pdiff(i)-cerror(i)],  'color', clr2, 'linewidth', lw); 
        h2 = plot(b1, [pdiff(i)+ferror(i) pdiff(i)-ferror(i)], 'color', clr1, 'linewidth', lw); 
    end
    
    plot(bls, pdiff, '.k', 'markersize', 20); 
    
    legend([h1, h2], 'cleared', 'forested', 'location', 'southwest'); 
    ylim([-4 4]); 
    xlabel('Perp Baseline (m)'); 
    ylabel('Phase Difference (radians)'); 
    title(['\Theta_{diff} Between Forested and Cleared Areas (' pol ')']); 
    box on
    f2 = gcf; 
    
    %print(f2, [datafol 'analysis/phsdiff_area' num2str(a) '_.pdf'], '-dpdf', '-painters');
    
% print(f1, [datafol 'analysis/bl_phs_' pol '_v' num2str(a) '.pdf'], '-dpdf', '-painters');
% print(f2, [datafol 'analysis/bl_phs_error_' pol '_v' num2str(a) '.pdf'], '-dpdf', '-painters');
% print(f1, [datafol 'analysis/bl_phs_' pol '_v' num2str(a) '.jpg'], '-djpeg');
% print(f2, [datafol 'analysis/bl_phs_error_' pol '_v' num2str(a) '.jpg'], '-djpeg');
    
    
    
    
    
    
    
%     
% close 
%           figure; 
%             subplot(1,2,1); 
%             histogram(p1, [-3.2:.1:3.2]); hold on; %i=11
%             xlim([-3.3 3.3]);
%             plot([p1m p1m], ylim, 'k-', 'linewidth', 2); 
%             plot([mean(p1) mean(p1)], ylim, 'k--', 'linewidth', 2); 
%             xlabel('phase values in box 1'); 
%             ylabel('number of pixels'); 
%             plot([p1m-e1 p1m-e1], ylim, 'r-'); 
%             plot([p1m+e1 p1m+e1], ylim, 'r-'); 
% 
%             subplot(1,2,2); 
%             histogram(p2, [-3.2:.1:3.2]); hold on; 
%             xlim([-3.2 3.2]);
%             plot([p2m p2m], ylim, 'k-', 'linewidth', 2); 
%             plot([mean(p2) mean(p2)], ylim, 'k--', 'linewidth', 2);
%             xlabel('phase values in box 2'); 
%             ylabel('number of pixels'); 
            %plot([p2m-e2 p2m-e2], ylim, 'k-'); 
            %plot([p2m+e2 p2m+e2], ylim, 'k-'); 
            
            
            
%             keyboard
%             % plot data 
%             close
%             figure; hold on; 
%             pcolor(phs); shading flat; hold on; 
%             plot([x1 x1 x2 x2 x1], [y1 y2 y2 y1 y1], 'k-');  
%             plot([xx1 xx1 xx2 xx2 xx1], [yy1 yy2 yy2 yy1 yy1], 'k-');  
%             axis([230 390 200 350])
%             colormap jet
%             linkaxes
    
    
    
    
    