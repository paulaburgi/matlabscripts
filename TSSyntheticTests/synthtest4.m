clear; 
n        = 30;    %number of interferograms
dz       = 3;      %height err/unit baseline
vel      = 10;      %velocity
numnoise = 2000;   %noise, each with own baseline
bpr      = 1;      %baseline range

ts       = linspace(0,1,n+1)';  %time vector
dt       = diff(ts);
tp       = ts(1:end-1)+dt/2; %t for plotting ints.
Gdt      = diag(dt); %make now, makes easier to remove values later.
nc       = floor(n/3);
ft       = [zeros(n,1);dz]; %time variable true model

tint     = and(ts>0.3,ts<0.65);
ft(tint) = vel;

%make bps
bps      = 2*bpr*rand(n+1,numnoise)-bpr;
bps      = diff(bps);

%impose smoothness on velocities
smoo     = -diag(ones(1,n-1),1)-diag(ones(1,n-1),-1);
ssum     = sum(smoo,2);
smoo     = smoo-diag(ssum);
%smoo=smoo(2:end-1,:);

lambda   = 0.05;

for i=1:n %loop over clearing times
    use    = [1:i-1 i+1:n]'; %don't use co-clearing ints.
 
    for j=1:numnoise
        Gt             = [Gdt(use,:) bps(use,j)]; %time-variable
        Gtc            = Gt;%greens function for clearing
        Gtc(use>i,end) = 0; %dates after clearing at ti have zero dem error
            
        %allow variable velocity
        Gtr            = [Gt;lambda*smoo zeros(size(smoo,1),1)]; %regularize
        Gg2            = inv(Gtr'*Gtr)*Gt';
        R2             = Gg2*Gtc;
        
        allm2(i,j,:)   = R2*ft;

    end
end


figure
is=round(n*[0.2 0.4 0.6 0.8]);
for i=1:length(is)
    m2std  = squeeze(std(allm2(is(i),:,:),[],2))';
    m2mean = squeeze(mean(allm2(is(i),:,:),2))';
    subplot(2,2,i)
    plot(tp,m2mean(1:n),'b');
    hold on
    plot([tp tp],m2mean(1:n)'*[1 1]+m2std(1:n)'*[-1 1],'b--')
    plot(tp,ft(1:n),'k')
    ax=axis;
    plot(tp(is(i))*[1 1],ax(3:4),'k:')
    xlabel('time')
    ylabel('inferred velocity')
    title(['logging at t=' num2str(tp(is(i)))])
end




stdall = 100*squeeze(std(allm2(:,:,1:end-1),[],2));
figure; 
imagesc(stdall); 
