% circ_phs_error
% Demonstrate relationship btw correlation value (R) and standard deviation

clear 
close all

pf_fol  = '/data/pmb229/isce/p222f870/';
datafol = [pf_fol 'data/']; 

% Define parameters
    ntrials  = 1000; 
    nsamples = 20*10;   % e.g. boxsize = 10x10, nsamples = 100  
    stdi      = 0:0.1:3; 

% Loop through std values
    R     = [];
    r_all = [];
    r_std = [];
    for i = 1:length(stdi)
        % define ntrials x nsamples matrix, for given std
            m = stdi(i).*randn(ntrials, nsamples); 

        % calculate mean and correlation (r) value using circular stats
            sum_sinm = sum(sin(m),2)./size(m,2); 
            sum_cosm = sum(cos(m),2)./size(m,2); 
            r     = sqrt((sum_sinm.^2) + (sum_cosm.^2)); 
            mm1  = sum_sinm./r; 
            mm2  = sum_cosm./r; 
            mm  = atan(mm1./mm2);

        % collect values for given stds
            R = [R; mean(r)]; 
            r_all = [r_all r]; 
            r_std = [r_std; std(r)];
    end

% equation for R vs std
    R_eq = 0:.01:1; 
    std_eq = sqrt(-2.*log(R_eq)); % source?? 

% plot results
    figure; hold on; box on; 
    for j=1:length(stdi)
        %plot(stdi(j), r_all(:,j), 'r.'); 
        %ya=minmax(r_all(:,j)); 
        %h1 = plot([stdi(j) stdi(j)], ya, '-', 'color', [.8 .8 .8]); 
    end
    h0 = plot(stdi, R+r_std, '--', 'color', lines(1)); 
    h0 = plot(stdi, R-r_std, '--', 'color', lines(1)); 
    h2 = plot(std_eq, R_eq, 'k', 'linewidth', 2);
    h3 = plot(stdi, R, 'color', lines(1), 'linewidth', 2); 

% label plot
    xlabel('Standard Deviation (\sigma)'); 
    ylabel('Correlation (R)'); 
    legend([h2, h3, h0], 'equation: \sigma = {\surd}( -2ln(R) )', ...
                         'for ntrials, mean of R for given std', ...
                         'for ntrials, std of R for given std', ... %'for ntrials, min to max of R values',...
                         'location', 'northeast'); 
    %set(gca, 'fontsize', 14); 
    axis([0 3.1 0 1]);
    title('Std vs R for box size 10 x 20 pixels'); 
    
    %print(gcf, [datafol 'analysis/gm_phserror.jpg'], '-djpeg');


