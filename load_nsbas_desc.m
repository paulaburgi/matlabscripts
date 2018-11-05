clear all
clc 
%close all

h5file='Stack/NSBAS-PARAMS.h5'; %mejor coherencia, 64 ints con 2/3 de los pixeles coherente
%h5file='NSBAS-PARAMS_2012_48.h5';
%h5file='NSBAS-PARAMS_2012_34.h5';
h5disp(h5file);
numints0=length(h5read(h5file,'/bperp'));

%change this as required
date0=datenum(2007,07,13);
%date0=datenum(2009,07,18);
intthresh=2/3; %this works just fine
%intthresh=0.5; %this works just fine
intthresh=round(numints0*intthresh); 

%dont change anything below this line
a=h5read(h5file,'/ifgcnt');
tims=double(h5read(h5file,'/tims'));
rawts=double(h5read(h5file,'/rawts'));
ifgcnt=double(h5read(h5file,'/ifgcnt'));

numints=size(rawts);
numints=numints(3);

%convert ints to vectors, so the linear approximation can be carried out
ints=[];
for i=1:numints
    temp=rawts(:,:,i);
    ints=[ints temp(:)];
end

ifgcnt=ifgcnt(:);
tic
for i=1:length(ints(:,1)) %iterate across all pixels
    vect=ints(i,:)/10; %mm to cm
    if ifgcnt(i)>=intthresh %number of ints larger than threshold
    %if nansum(vect)>intthresh
        k = ~isnan(vect); %get rid of nans
%        p = polyfit(tims(k),vect(k)',1); %fit a pol to non NaNs
        p = [tims(k) ones(size(tims(k)))]\vect(k)';
        ratemap(i) = p(1);
        linearcoef(i) = p(2);
        residual = vect(k)' - [tims(k) ones(size(tims(k)))]*p;
        res(i) = sqrt(residual'*residual/length(residual));
    else
        ratemap(i)=NaN;
        linearcoef(i) = NaN;
        res(i) = NaN;
    end
end
toc

ratemap=reshape(ratemap,[size(rawts(:,:,1))]);
ratemap=flipud(ratemap');
res=reshape(res,[size(rawts(:,:,1))]);
res=flipud(res');
nanvar(res(:))
linearcoef=reshape(linearcoef,[size(rawts(:,:,1))]);
linearcoef=flipud(linearcoef');

for i=1:numints
    Rts(:,:,i)=flipud(rawts(:,:,i)');
end

maxvel=find(ratemap==(max(max(ratemap))));
[yp,xp] = ind2sub(size(ratemap),maxvel); %yp = row, xp = column
%  yp=1078;
%  xp=2070;
% xp=471; 
% yp=30; 
% xp = 2053; % 0
% yp = 1889; 
% xp = 1822; 
% yp = 861; % +
xp = 1837; 
yp = 877; % -
[xp yp]

%mask errors
%ratemap(655:685,704:735)=NaN;

%% plot
figure
    hold on
    pcolor(ratemap); shading flat;
    plot(xp,yp,'.k','MarkerSize',20)
    axis equal
    %axis([300 900 350 950])
    colorbar
    caxis([0 max(ratemap(:))])
    caxis([-30 30])
    colormap jet

wvl=5.54657;
wvl=10;
ratemap_wrap=wrapToPi((ratemap-min(ratemap(:)))*4*pi/wvl);
ratemap_wrap=wrapToPi((ratemap)*4*pi/wvl);

% figure
%     hold on
%     pcolor(ratemap_wrap)
%     plot(xp,yp,'.k','MarkerSize',20)
%     shading flat
%     axis equal
%     %axis([300 900 350 950])
%     colorbar
%     title([num2str(numints0), ' interferograms, ',num2str(intthresh/numints0*100) ' of coherent ints'])

ratemap_wrap=(ratemap_wrap/(2*pi) +1/2 ) * wvl/2; %wrap to wvl/2
ratemap_wrap(isnan(ratemap_wrap))=-99;

for i=1:numints
    maxdef(i)=[Rts(yp,xp,i)]/10;
end
% for i=1:numints
%     maxdef2(i)=[Rts(450,1000,i)]/10;
% end
% lfit2=tims*ratemap(450,1000) + linearcoef(450,1000)

lfit=tims*ratemap(yp,xp) + linearcoef(yp,xp);
rms_lin=sqrt(dot(lfit-maxdef',lfit-maxdef')/length(lfit-maxdef') );

tims=tims*365.24+date0;

figure
    hold on
    box on
    plot(tims,maxdef,'o','MarkerSize',8,'MarkerFaceColor',[0 0.15 0.5],'MarkerEdgeColor','k')
    plot(tims,lfit,'LineWidth',2,'Color',[1 1 1]*0.15)
%     plot(tims,maxdef2,'o','MarkerSize',8,'MarkerFaceColor',[1 0.25 0],'MarkerEdgeColor','k')
%     plot(tims,lfit2,'LineWidth',2,'Color',[1 1 1]*0)
    xlabel('Time (years)')
    ylabel('LOS disp (cm)')
    datetick('x','yyyy-mm')
    %print -depsc2 ts_perf.eps
    save max_uplift_desc tims maxdef
    
[nx,ny]=size(rawts(:,:,1));
delta=1; %displace CSK by 1 row to the S (error when manually merging SRTM frames), so it will match RS2
rr=ratemap(1+delta:ny,1:nx);
ratemap(isfinite(ratemap))=NaN;
ratemap(1:ny-delta,1:nx)=rr;
rr2=ratemap_wrap(1+delta:ny,1:nx); 
ratemap_wrap(isfinite(ratemap_wrap))=NaN;
ratemap_wrap(1:ny-delta,1:nx)=rr2;
yp=yp-1;
    
sAs = single([flipud(0*ratemap) flipud(ratemap/100)]'); %save grid in m
fid = fopen('csk_desc.geo','wb');
fwrite(fid,sAs,'real*4');
fclose(fid);


ratemap2=ratemap;
ratemap2(isnan(ratemap2)==1)=-99;

%ginput coords for plot in GMT
nx=size(ratemap,2);
ny=size(ratemap,1);

x1=-123.76027777777779;
y1=43.5627777777777799346;
dx=0.0002777777777777778;

y2=y1+dx*(ny-1);
x2=x1+dx*(nx-1);
x=x1:dx:x2;
y=flipud(y1:dx:y2);
disp([num2str([x1+dx*(xp-1) y1+dx*(yp-1)])])
max_uplift_coords=[x1+dx*(xp-1) y1+dx*(yp-1)];
%save max_uplift_coords_desc.txt max_uplift_coords -ascii
grdwrite2(x,y,ratemap2,'test4.grd')

% 
% ind2crop=52; %crop to get rid of the NaN's beneath the DEM
% ratemap2=ratemap2(ind2crop:end,:);
% y=y(ind2crop:end);
% disp([min(x) max(x) min(y) max(y)])
% 
% temp = fliplr(ratemap2');
% temp1 = temp(:);
% fid  = fopen('csk_desc.r4','w');
% fwrite(fid, temp1, 'real*4');
% fclose(fid);