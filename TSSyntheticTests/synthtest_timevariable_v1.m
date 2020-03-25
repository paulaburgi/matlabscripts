%keyboard
clear

sr       = 9.48e5;   % Slant Range (m)
los      = 43.0;     % Line-of-sight angle (degrees)
l        = 0.055;    % Wavelength (m)
r2m      = (4*pi)/l; % Radians to meters conversion
n        = 10;       % number of interferograms
dz       = 30;       % height err/unit baseline
vel      = 0.05;     % velocity (m/year)
numnoise = 100;    % noise, each with own baseline
bpr      = 60;       % baseline range
t        = 12;       % repeat interval (days)


ts       = (0:t:n*t)';         % time vector
dt       = diff(ts);
tp       = ts(1:end-1)+dt/2;   % t for plotting ints.
Gdt      = r2m*diag(dt)/365;   % make now, makes easier to remove values later.
ft       = [zeros(n,1);dz];    % time variable true model

tint     = and(ts>ts(end)*0.3,ts<ts(end)*0.65);
%tint     = and(ts>ts(end)*0.2,ts<ts(end)*0.4);
ft(tint) = vel;

%make bps
bps      = 2*bpr*randn(n+1,numnoise)-bpr;
bps      = diff(bps);
bps      = (r2m*bps)./(sr*sind(los));

%impose smoothness on velocities
smoo     = -diag(ones(1,n-1),1)-diag(ones(1,n-1),-1);
ssum     = sum(smoo,2);
smoo     = smoo-diag(ssum);
%smoo=smoo(2:end-1,:);
smoo(1,:)   = [2 -1 zeros(1,size(smoo,1)-3) -1]; 
smoo(end,:) = [-1 zeros(1,size(smoo,1)-3) -1 2]; 
smoo = smoo*-1; 

allm      = {};
lambdai   = [1.4,5]; %1e-5;% 1.4; 4;
% lambdai   = logspace(-5,2, 20); 
ttt = 0; 
ddd = length(lambdai)*n; 

for k=1:length(lambdai)
for i=1:n                    % Loop over clearing times
    use    = [1:i-1 i+1:n]'; % Don't use co-clearing ints.
 
    for j=1:numnoise
        Gt             = [Gdt(use,:) bps(use,j)]; % Time-variable
        Gtc            = Gt;                      % Greens function for clearing
        Gtc(use>i,end) = 0;                       % Dates after clearing at ti have zero dem error
            
        % Allow variable velocity
        lambda         = lambdai(k); 
        L              = [smoo zeros(size(smoo,1),1)]; 
        Gtr            = [Gt;lambda*L];    % Regularize
        Gg2            = inv(Gtr'*Gtr)*Gt';
        R2             = Gg2*Gtc;
        
        m              = R2*ft;
        allm2(i,j,:)   = m(1:end-1);
        d              = Gtc*ft; 
        rn(i,j)        = norm((Gt*m - d)); % Residual norm
        mn(i,j)        = norm(L*m);        % Model norm
    end
    ttt = ttt + 1; 
    nnn = datestr(now); 
    disp([num2str(ttt) '/' num2str(ddd) ', ' nnn(13:end)]); 
end
    rna(k) = mean(mean(rn));
    mna(k) = mean(mean(mn)); 
    allm   = [allm; allm2];
    % disp(num2str(k)); 
end
%%
% figure;
% plot(log10(rna(1:5)), log10(mna(1:5)), '.-')

close all; 
figure('units', 'normalized', 'outerposition', [.3 .8 .3 1]); 
is=round(n*[0.2 0.4 0.6 0.8]);
% is=round(n*[0.02:0.03:.3]);
% is=round(n*[0.1:0.1:1]);
loc = [0.7 0.5 0.3 0.1]; 
lg  = ['A'; 'B'; 'C'; 'D']; 

idx = find(tint~=0); 
tintf = [tint(1:idx(1)-1); tint(idx(1)); tint(idx); tint(idx(end)); tint(idx(end)+1:end)]*vel; 
tpf   = [tp(1:idx(1)-1); tp(idx(1))-(dt(1)/2); tp(idx); tp(idx(end))+(dt(1)/2); tp(idx(end)+1:end)]; 

nplt = length(is); 
cc=[     0    0.4470    0.7410; ...
    0.8500    0.3250    0.0980; ...
    0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560]; 

for i=1:length(is)
    allm22 = cell2mat(allm(2));     
    m2std2  = 100*squeeze(std(allm22(is(i),:,:),[],2))';
    m2mean2 = 100*squeeze(mean(allm22(is(i),:,:),2))';
    allm2 = cell2mat(allm(1)); 
    m2std  = 100*squeeze(std(allm2(is(i),:,:),[],2))';
    m2mean = 100*squeeze(mean(allm2(is(i),:,:),2))';
%     subplot(2,nplt/2,i); hold on; box on; 
    %subplot(4,1,i); hold on; box on; 
    subplot('Position', [0.1 loc(i) 0.8 0.2]); 
    set(gca, 'FontName', 'Arial', 'fontsize', 14); hold on;  hold on; box on; 
    plot(tp(is(i))*[1 1],[-23 23],'--k'); 
    %plot(tp,100*ft(1:n),'color', [0.7 0.7 0.7], 'linewidth', 2)
    stairs(tpf,100*tintf(1:n+2), '.-', 'color', [0.53 0.53 0.53], 'linewidth', 2.5)
    %bar(tp, 100*ft(1:n), 'edgecolor', 'none', 'facecolor', [0.6 0.6 0.6]); 
    plot(tp,m2mean(1:n), '.-','color',  cc(1,:), 'linewidth', 1.5); 
    plot(tp,m2mean2(1:n), '.-','color',  cc(4,:), 'linewidth', 1.5); 
    fill([(tp); flipud(tp)], [(m2mean(1:n)+m2std(1:n)),fliplr(m2mean(1:n)-m2std(1:n))], cc(1,:), 'edgecolor', 'none', 'LineStyle', '-', 'FaceAlpha', 0.2);
    fill([(tp); flipud(tp)], [(m2mean2(1:n)+m2std2(1:n)),fliplr(m2mean2(1:n)-m2std2(1:n))], cc(4,:), 'edgecolor', 'none', 'LineStyle', '-', 'FaceAlpha', 0.2);
    xlim([ts(1)-12 ts(end)+12]); 
    %ylim([ax(3) ax(4)]);
    ylim([-23 23]);
    if i == 1
        ax=axis;
        %title(['lambda: ' num2str(lambda)]); 
    end
    if i == 4
        xlabel('Time (days)')
    elseif i == 2
        ylabel('Inferred Velocity (cm/yr)                                 '); 
        set(gca, 'XTicklabel', []); 
    else
        set(gca, 'XTicklabel', []); 
    end
    
    yticks([-15 0 15]); 
    %text(335, 19, lg(i), 'fontsize', 15, 'fontname', 'arial'); 
    
    grid on; 
    xlim([tp(1) tp(end)]); 
    %title({['Logging on day: ' num2str(tp(is(i)))]})
    
    
end
linkaxes; 
% keyboard


% stdall = 100*squeeze(std(allm2(:,:,:),[],2));
% mean(mean(stdall))
% figure; 
% imagesc(stdall); 
% title(['lamda: ' num2str(lambda)])











