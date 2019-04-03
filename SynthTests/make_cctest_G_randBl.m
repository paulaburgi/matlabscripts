function [G, zrdi] = make_cctest_G_randBl(dc, du)
    dta    = diff(du); 
    nvels  = length(du)-1;
    nints  = length(dc); 

% build G matrix
    G  = zeros(nints, nvels);
    for i = 1:size(G, 1)
        x   = find(du == dc(i,1)); 
        y   = find(du == dc(i,2)); 
        
        dtj = dta(x:y-1);
        for j = 1:length(dtj)
            G(i,x+j-1) = dtj(j); 
        end
    end
    d2y = 0.00274; % convert days to years
    Ga  = G.*d2y; 
    
% Add rows s.t. mean vel = avg vel for each unconstrained date
    rdi    = find(sum(G)' == 0);
    nrdi   = length(rdi);
    avgvel = (1/(nvels-1)); %-nrdi)); 
    G      = [Ga; ones(nrdi, nvels)*avgvel]; 
    rng    = 1:nrdi; 
    for i = rng
        idx = i-1; 
        G(end-idx, rdi(i)) = -1; 
    end
    zrdi = zeros(nrdi, 1); % need to add to end of other vectors

 % get rid of 1st column s.t. 1st time step = 0
    G    = G(:,2:end); 

    
    
    
    
    
    
    
    
    
    
    
    