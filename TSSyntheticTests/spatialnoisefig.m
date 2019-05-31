nx       = 100;
x        = linspace(0,1,nx);
ns       = 0.5; %number of sinusoids for "ramp" - small = more linear
Var      = 0.1; %noise level for random white noise

nd       = 10; %number of dates
postd    = nd*0.5; %timing of clearcut
poststat = (1:nd)>=postd;

bp       = linspace(-1,1,nd);
bp       = bp(randperm(nd)); %random, but evenly distributed baselines

%random values for ramp
sinphase = rand(1,nd)*2*pi;
sinmag   = 0.75*randn(1,nd);

%mark out boxcars;
%nb=3;sigbox=ceil(x*nb)-round(x*nb); ;uniform boxcars
b        = [0.1 0.2 ;0.4 0.55; 0.75 0.9];
bcar     = zeros(1,nx);
for i=1:3
    bcar(and(x>b(i,1),x<=b(i,2)))=0.25+i/4; %slightly different sizes for clarity in last figure
end

%after time postd, remove 3rd boxcar
bcarpost=bcar;
bcarpost(and(x>b(3,1),x<=b(3,2)))=0;  


%points: perm clearcut, middle tree boxcar, third boxcar that gets cleared
xp=[70 50 80];



dx   = 10; %spatial scale for intra-clearcut filter; note, this is in pixels, not distance
kern = exp(-(-dx:dx).^2/2/(dx/3)^2);

for i=1:nd
    allnoise(i,:)  = Var*randn(1,nx);
    allramp(i,:)   = sinmag(i)*sin(x*2*pi*ns+sinphase(i));
    
    %signal from boxcars, pre or post
    allsigbox(i,:) = bcar*bp(i);
    if(i<postd)
        allsigbox2(i,:)=allsigbox(i,:);
    else
        allsigbox2(i,:)=bcarpost*bp(i);
    end
    
    %if clearcuts known (sigbox=0), filter out atm?
    tmp        = allsigbox2(i,:)+allnoise(i,:)+allramp(i,:); % signal with random noise and ramp
    good       = bcar==0;
    tmp(~good) = 0;
    a          = conv(tmp,kern,'same');
    b          = conv(good,kern,'same');
    allsigbox2_filt(i,:)=a./b;
end
allsignoise=allsigbox2+allnoise;
allsignoiseramp=allsignoise+allramp;


%nice colorscale from red to blue
nd2=ceil(nd/2);
colors=[ones(nd2,1) linspace(0,0.9,nd2)' linspace(0,0.9,nd2)';linspace(0.9,0,nd2)' linspace(0.9,0,nd2)' ones(nd2,1);];


close all; 
figure
subplot(2,3,1)
for i=1:nd
    plot(x,allsigbox(i,:),'-','color',colors(i,:),'linewidth',abs(bp(i))*3)
    hold on
end
plot(x(xp),0,'ko','markerfacecolor','k')
text(x(xp),0*xp,num2str([1:3]'))
axis([0 1 -2 2])
title('forest+clearcuts, random baselines')
xlabel('distance'); 
ylabel('signal'); 


subplot(2,3,2)
for i=1:nd
    plot(x,allsigbox2(i,:),'-','color',colors(i,:),'linewidth',abs(bp(i))*3)
    hold on
end
axis([0 1 -2 2])
title('new clearcut halfway')


subplot(2,3,3)
for i=1:nd
    plot(x,allsignoise(i,:),'-','color',colors(i,:),'linewidth',abs(bp(i))*3)
    hold on
end
axis([0 1 -2 2])
title('add random noise')

subplot(2,3,4)
for i=1:nd
    plot(x,allsignoiseramp(i,:),'-','color',colors(i,:),'linewidth',abs(bp(i))*3)
    hold on
end
axis([0 1 -2 2])
title('add longer wavelength atm')

subplot(2,3,5)
for i=1:nd
    plot(x,allsignoiseramp(i,:)-allsigbox2_filt(i,:),'-','color',colors(i,:),'linewidth',abs(bp(i))*3)
    hold on
end
axis([0 1 -2 2])
title('removed filtered clearcut signal from total')

subplot(2,3,6)
for i=1:nd
    plot(x,allramp(i,:)-allsigbox2_filt(i,:),'-','color',colors(i,:),'linewidth',abs(bp(i))*3)
    hold on
end
axis([0 1 -2 2])
title('diff between panels c and e')




%%%now plot baseline vs. phs
ls=[3 2 1]; %so plots don't overla
[jnk,bpid]=sort(bp);
figure
subplot(2,3,1)
hold on
for i=1:3
    plot(bp(bpid),allsigbox(bpid,xp(i)),'-','linewidth',ls(i))
end
plot(bp(poststat),allsigbox(poststat,xp(3)),'k.') %plot points in third clearcut at times post-clearcutting
title('constant clearcuts');

subplot(2,3,2)
hold on
for i=1:2
    plot(bp(bpid),allsigbox2(bpid,xp(i)),'-','linewidth',ls(i))
end
plot(bp(~poststat),allsigbox2(~poststat,xp(3)),'.','markersize',20) %plot points in third clearcut at times pre-clearcutting
plot(bp(poststat),allsigbox2(poststat,xp(3)),'k.') %plot points in third clearcut at times post-clearcutting
title('time-variable clearcut at pt 3')
subplot(2,3,3)
hold on
for i=1:2
    plot(bp(bpid),allsignoise(bpid,xp(i)),'-','linewidth',ls(i))
end
plot(bp(~poststat),allsignoise(~poststat,xp(3)),'.','markersize',20) %plot points in third clearcut at times pre-clearcutting
plot(bp(poststat),allsignoise(poststat,xp(3)),'k.') %plot points in third clearcut at times post-clearcutting

title('plus noise')

subplot(2,3,4)
hold on
for i=1:2
    plot(bp(bpid),allsignoiseramp(bpid,xp(i)),'-','linewidth',ls(i))
end
plot(bp(~poststat),allsignoiseramp(~poststat,xp(3)),'.','markersize',20) %plot points in third clearcut at times pre-clearcutting
plot(bp(poststat),allsignoiseramp(poststat,xp(3)),'k.') %plot points in third clearcut at times post-clearcutting

title('noise+ramp')
subplot(2,3,5)
hold on
for i=1:2
    plot(bp(bpid),allsignoiseramp(bpid,xp(i))-allsigbox2_filt(bpid,xp(i)),'-','linewidth',ls(i))
end
plot(bp(~poststat),allsignoiseramp(~poststat,xp(3))-allsigbox2_filt(~poststat,xp(3)),'.','markersize',20) %plot points in third clearcut at times post-clearcutting
plot(bp(poststat),allsignoiseramp(poststat,xp(3))-allsigbox2_filt(poststat,xp(3)),'k.') %plot points in third clearcut at times post-clearcutting

title('rem filt clearcut from total')




close all; 
% figure for andes
figure
subplot(2,2,1)
for i=1:nd
    plot(x,allsignoise(i,:),'-','color',colors(i,:),'linewidth',abs(bp(i))*3)
    hold on
end
axis([0 1 -2 2])
title('Data with random noise')
xlabel('distance'); 
ylabel('signal'); 

subplot(2,2,2)
for i=1:nd
    plot(x,allsignoiseramp(i,:),'-','color',colors(i,:),'linewidth',abs(bp(i))*3)
    hold on
end
axis([0 1 -2 2])
title('Data with random noise and ramp')
xlabel('distance'); 
ylabel('signal'); 

subplot(2,2,3)
for i=1:nd
    plot(x,allsignoiseramp(i,:)-allsigbox2_filt(i,:),'-','color',colors(i,:),'linewidth',abs(bp(i))*3)
    hold on
end
axis([0 1 -2 2])
title('HPF data with random noise and ramp')
xlabel('distance'); 
ylabel('signal'); 

subplot(2,2,4)
for i=1:nd
    plot(x,allramp(i,:)-allsigbox2_filt(i,:),'-','color',colors(i,:),'linewidth',abs(bp(i))*3)
    hold on
end
axis([0 1 -2 2])
title('Diff btw a & c')
xlabel('distance'); 
ylabel('signal'); 