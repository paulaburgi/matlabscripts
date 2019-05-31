dn_all = [1:16]'; 
gidx   = [[1:15]' [2:16]'; 4 6; 6 8; 10 12];
dn     = dn_all(gidx); 
    
% bli   = [1.4969; 1.0182; 1.2191; 1.5412; 0.5108; 1.0643; 0.6119; 0.9239; 1.3436; ...
%          0.3724; 1.4255; 1.2053; 0.4849; 1.6157; 0.1; 1.0970; 0.7212; 0.8430];
bli   = [0; 0.51; -0.53; 0.24; -0.25; 0.05; 0.62; 0.45; -0.63; -0.23; 0.33; -0.28; -0.39; 0.14; 0.55; -0.32];
% bli   = [0; (rand(length(dn_all),1))-0.5]; 
blii  = bli(gidx); 
blia  = blii(:,2)-blii(:,1); 
%bli   = blia+abs(min(blia)); 
bli = abs(blia); 
G     = bli;

ints  = bli; 

% test 1
ints1 = zeros(length(ints),1); 
m1    = inv(G'*G)*G'*ints1; 

% test 2
ints2 = [ints(1:3); zeros(13,1)];
G2    = [G(1:3); G(5:15); G(17:end)]; 
m2    = inv(G2'*G2)*G2'*ints2; 

% test 3
ints3 = [ints(1:10); zeros(4,1); ints(end-2:end-1)];
G3    = [G(1:10); G(12:15); G(end-2:end-1)]; 
m3    = inv(G3'*G3)*G3'*ints3; 

% test 4
ints4 = ints; 
G4    = G; 
m4    = inv(G4'*G4)*G4'*ints4; 

close all; 
figure('units', 'normalized', 'outerposition', [.1 .6 .2 .9]); hold on; box on;
subplot(4,1,1); hold on; box on; 
plot(G, ints1, '.', 'markersize', 15, 'color', [0.6 0 0])
plot([0; G], [0; G*m1], 'color', [0.6 0 0]);
ylim([-0.040 max(bli)+.21]); 
subplot(4,1,2); hold on; box on; 
plot(G2, ints2, '.', 'markersize', 15, 'color', [0.6 0 0])
plot([0; G2], [0; G2*m2]);
ylim([-0.040 max(bli)+.21]); 
subplot(4,1,3); hold on; box on; 
plot(G3, ints3, '.', 'markersize', 15, 'color', [0.6 0 0])
plot([0; G3], [0; G3*m3]);
ylim([-0.040 max(bli)+.21]); 
subplot(4,1,4); hold on; box on; 
plot(G4, ints4, '.', 'markersize', 15, 'color', [0.6 0 0])
plot([0; G4], [0; G4*m4]);
ylim([-0.040 max(bli)+.21]); 

figure('units', 'normalized', 'outerposition', [.3 .6 .65 .3]);; hold on; box on; 
plot(dn, blii, '.', 'markersize', 15, 'color', [0.4 0.4 0.4]); 
for i = 1:length(ints)
    plot(dn(i,:), blii(i,:), 'color', [0.4 0.4 0.4]); 
end
xlim([0 dn_all(end)+1]); 
ylim([min(min(blii))-.21 max(max(blii))+.21]); 


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




