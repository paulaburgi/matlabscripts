% see if I can reproduce trend in cctest_avgvel_randBl.m
clear
close all; 

% slant range, line of sight, dem error
sr         = 9.48e5; 
los        = 43.0; 
l          = 0.055; 
% l           = 0.236; 
% sr          = 8.5e5; 
% los         = 38.7;

dz         = 30; 
lp = l/(4*pi); 
pl = (4*pi)/l; 


Vavg_all = []; 
for j = 1:1000
    
    nints = 30; 
    ints  = zeros(nints, 1);
    
    Bl1    = [0 randn(1,1); randn(nints-1,2)]; 
    Bl1    = diff(Bl1')'; 
    
    Bl2    = [0; randn(nints,1)];
    Bl2    = Bl2([[1:nints]' [2:nints+1]']);
    Bl2    = [diff(Bl2')']; 
    
    
    
    
    Bl    = Bl1; 
    Bl0    = Bl*70*100; %(Bl * m->cm * Bp scale) % sentinel:70, alos:300
    %Bl    = (2*Bl0)./(sr.*sind(los)); 
    %Bl    = (l*Bl0)./(4*pi*sr*sind(los)); 
    Bl    = (4*pi*Bl0)./(l*sr*sind(los)); 
    G     = [ones(nints, 1)*0.00274*12 Bl];

    Vavg = [];
    for i = 0:nints
        if i == 0
            Gi    = G; 
            intsi = ints; 
        elseif i == nints
            ints(i) = Bl(i)*dz;
            intsi   = ints; 
            Gi = G; 
        elseif i ~= 0 && i ~= nints
            ints(i) = Bl(i)*dz;
            intsi   = ints([1:i-1 i+1:end]); 
            Gi      = G([1:i-1 i+1:end], :); 
        end
        Ginv = inv(Gi'*Gi)*Gi'; 
        m    = Ginv*intsi; 
        Vavg = [Vavg; m(1)]; 
    end

    Vavg_all = [Vavg_all Vavg]; 
end

close all; 
figure; hold on; box on; 
plot(0:nints+4, zeros(nints+5,1), 'color', [0.7 0.7 0.7], 'linewidth', 2); 
Vm = mean(Vavg_all,2);
Vs = std(Vavg_all')';
errorbar(Vm, Vs, 'k'); 













































