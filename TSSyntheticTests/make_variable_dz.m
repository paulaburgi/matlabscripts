function dz = make_variable_dz(dz0, dc, gr)
    % dz0 = 30; 
    % gr  = 0.00274*2; % growth rate of trees per day (0.00274 m/day = 1 m/year)

t1  = min(min(dc)); 
t2  = max(max(dc)); 

t_all  = t1:t2; 
ddz    = ones(length(t_all),1).*gr;  
dz_all = cumsum([dz0; ddz]); 

dz = zeros(length(dc),1); 
for i = 1:length(dc)
    idx1  = find(dc(i,1) == t_all); 
    idx2  = find(dc(i,2) == t_all); 
    dz(i) = mean([dz_all(idx1) dz_all(idx2)]); 
end