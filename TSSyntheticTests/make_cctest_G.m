function [Gbl, zrdi, blg] = make_cctest_G(dc, du, bl, rparams)

    dta    = diff(du); 
    nvels  = length(du)-1;
    nints  = length(dc); 

    xa = [];
    ya = [];
% build G matrix
    G  = zeros(nints, nvels);
    xa = [];
    ya = [];
    for i = 1:size(G, 1)
        xi  = dc(i,1); 
        yi  = dc(i,2); 
        x   = find(du == xi); 
        y   = find(du == yi); 
        
        dtj = dta(x:y-1);
        for j = 1:length(dtj)
            G(i,x+j-1) = dtj(j); 
        end

        xa = [xa; x];
        ya = [ya; y];
    end
    d2y = 0.00274; % convert days to years
    Ga  = G.*d2y; 
    
% Add rows s.t. mean vel = avg vel for each unconstrained date
    rdi    = find(sum(G)' == 0);
%     if rdi ==1
%         rdi = [2];
%         Ga(1,:) = 0; 
%     end
    nrdi   = length(rdi);
    avgvel = (1/(nvels-1)); %-nrdi)); 
    G      = [Ga; ones(nrdi, nvels)*avgvel]; 
    rng    = 1:nrdi; 
    for i = rng
        idx = i-1; 
        G(end-idx, rdi(i)) = -1; 
    end
    zrdi = zeros(nrdi, 1); % need to add to end of future vectors

% add baseline term to G matrix
    %blg    = [(4.*pi.*bl)./(rparams(1).*rparams(2).*sind(rparams(3))); zrdi]; 
    blg    = [(2*bl)./(rparams(2).*sind(rparams(3))); zrdi]; 
    Gbl    = [G(:,2:end) blg]; 
    %Gbl    = [G blg]; 
    %Gbl    = [zeros(length(G),1) G blg]; 

    
    
    
    
    
    
    
    
    
    
    
    