parpool(10)
home=pwd;
dirs={'T130_T28_resamp'};
%dirs=dirs(1:2);
pol='_VV';

nt=length(dirs);
nd=0;
dates=[];
for i=1:length(dirs)
    tmp=dir([dirs{i} '/rel*.cor.geo']);
    %make a list of T28 dates
    for j=1:length(tmp)
        t=tmp(j).name;
		%chcek to see if date is a T28 date. Only do next two lines if so.
        dates(end+1).name=[tmp(j).folder '/' t];
        dates(end).dn=datenum(t(5:12),'yyyymmdd');
    end
end
nd=length(dates);
dn    = [dates.dn];
[jnk,sortid]=sort(dn);
dates = dates(sortid);
dn    = [dates.dn];



%T130 T28
nx_geo=822;
ny_geo=2538;


rdir = ['results_TS_T28_only' pol '/'];
rdate  = {'20170325','20170506','20170717','20171009','20171220','20180302','20180524','20181011','20181109'};
if(~exist(rdir,'dir'))
    mkdir(rdir)
end

dnr    = datenum(rdate,'yyyymmdd');
rdate  = rdate(dnr<max(dn));
dnr    = dnr(dnr<max(dn));
for i=1:length(dnr)
    rain(i).dnr=dnr(i);
    if(i<length(dnr))
        afid=find(and(dn>=dnr(i),dn<dnr(i+1)));
    else
        afid=find(dn>=dnr(i));
    end
    [jnk,sortid]=sort(dn(afid));
    rain(i).valid=afid(sortid);
end

for i=1:nd
    fidr(i)=fopen(dates(i).name,'r');
end

for i=1:length(dirs)
files=dir([dirs{i} '/*c0.cor.geo']);
for j=1:length(files)
    	tmp=[dirs{i} '/' files(j).name];
    fid0(j)=fopen(tmp,'r');
    end
end
%find last line written
testfile=[rdir rdate{1} '.mag0'];
if(exist(testfile,'file'))
    fid=fopen(testfile,'r');
    for j=1:ny_geo
        [tmp,count]=fread(fid,nx_geo,'real*4');
        if(count<nx_geo)
            break
        end
    end
    fclose(fid);
    
    online=j-1
    if(online<1)
        disp([ testfile  ' empty?'])
        return
    end
    for i=1:length(dnr)
        fidmag(i)  = fopen([rdir rdate{i} '.mag0'],'a');
        fidt(i)    = fopen([rdir rdate{i} '.time0'],'a');
        fidmagl(i)  = fopen([rdir rdate{i} '.maglow'],'a');
        fidmagh(i)  = fopen([rdir rdate{i} '.maghigh'],'a');
        fidte(i)    = fopen([rdir rdate{i} '.timeerr'],'a');
        fidn(i)     = fopen([rdir rdate{i} '.count'],'a');
 
    end
    
    fidres  = fopen([rdir 'resn0'],'a');
    
    for i=1:nd
        fseek(fidr(i),online*nx_geo*4,-1);
    end
    for i=1:length(fid0)
        fseek(fid0(i),online*nx_geo*4,-1);
    end
else
    online=0;
    for i=1:length(dnr)
        fidmag(i)  = fopen([rdir rdate{i} '.mag0'],'w');
        fidt(i)    = fopen([rdir rdate{i} '.time0'],'w');
        fidmagl(i)  = fopen([rdir rdate{i} '.maglow'],'w');
        fidmagh(i)  = fopen([rdir rdate{i} '.maghigh'],'w');
        fidte(i)    = fopen([rdir rdate{i} '.timeerr'],'w');
        fidn(i)     = fopen([rdir rdate{i} '.count'],'w');
    end
    fidres  = fopen([rdir 'resn0'],'w');
end

ft = fittype( 'a*exp(-b*x)', 'independent', 'x', 'dependent', 'y' );

%start
for j=online+1:ny_geo
    j
    dat=zeros(nd,nx_geo);
    for i=1:nd
        [tmp,count]=fread(fidr(i),nx_geo,'real*4');
        if(count>0)
            dat(i,1:count)=tmp;
        end
    end
    dat(dat==-9999)=NaN;
    dat(isinf(dat))=NaN;
    dat(dat==0)=NaN;
    
    count   = sum(isfinite(dat));
    
    c0s=nan(length(fid0),nx_geo);
    for i=1:length(fid0)
        [tmp]=fread(fid0(i),nx_geo,'real*4');
            c0s(i,:)=tmp;
    end
    c0s(isinf(c0s))=NaN;
    c0s(c0s==0)=NaN;
    c0s(c0s==-9999)=NaN;
    c0s=median(c0s,1);
    c=1-c0s;
    c0err=0.067*c-0.34*c.^2;
    
    
    goodid  = find(and(count>50,c0s>0.3));
   
    
    mags   = nan(length(dnr),length(goodid));
    times  = mags;
    maglow = mags;
    maghig = mags;
    timerr = mags;
    count  = mags;

    for k=1:length(dnr)
        valid  = [rain(k).valid];
        x=[dn(valid)-dnr(k)]';
        d=-log(dat(valid,:));  
        
        r=isfinite(d);
        tmp=[1:length(x)]'*ones(1,nx_geo).*r;
        tmp(tmp==0)=inf;
        ids=min(tmp);
        goodi=isfinite(ids);
        ids=sub2ind(size(d),ids(goodi),find(goodi));
        dfirst=zeros(1,nx_geo);
        dfirst(goodi)=d(ids);
        dfirst(~isfinite(dfirst))=0;
        
        
        dbefall=-log(dat(1:valid(1)-1,:)); %previous data points
        r=isfinite(dbefall);
        tmp=[1:valid(1)-1]'*ones(1,nx_geo).*r; %this makes a matrix of indices, with nans masked out
        befids=max(tmp);
        goodbef=befids>0;
        befids=sub2ind(size(dbefall),befids(goodbef),find(goodbef)); 
        dbef=zeros(1,nx_geo);
        dbef(goodbef)=dbefall(befids);%pick last non-nan before date, for each pixel in row
        dbef(~isfinite(dbef))=0; %just in case
        
        parfor i=1:length(goodid)
            myopt = fitoptions( 'Method', 'NonlinearLeastSquares','TolFun',1e-2);
            myopt.Display = 'Off';
            myopt.Lower = [0 0];
            
            data        = d(:,goodid(i));
            firstdif    = min(dfirst(goodid(i)),dfirst(goodid(i))-dbef(goodid(i)),'omitnan');
            gi          = find(isfinite(data));
            count(k,i)  = length(gi);
                         
                 
            if(length(gi)>=4)
                if(firstdif>c0err(goodid(i))*2)
                    c       = 1-dat(valid,goodid(i));
                    derr    = 0.067*c-0.34*c.^2;
                    dup     = c+derr;
                    ddn     = c-derr;
                    dupl    = -log(1-dup);
                    ddnl    = -log(1-ddn);
                    weights = sqrt(abs(1./(dupl-ddnl)/2));
                    weights(isinf(weights))=100;
                    
                    % Fit model to data.
                    myopt.StartPoint = [dfirst(goodid(i)) 0.1];
                    myopt.Weights    = weights(gi);
                    [fitresult]      = fit(x(gi),data(gi), ft, myopt );
                    results          = coeffvalues(fitresult);
                    confs            = confint(fitresult);
                    
                    if(sum(confs(:)<0)==0)
                        mags(k,i)   = exp(-results(1));
                        times(k,i)  = 1./results(2);
                        maglow(k,i) = exp(-confs(2,1));
                        maghig(k,i) = exp(-confs(1,1));
                        timerr(k,i) = -diff(1./confs(:,2));
                    else
                        maglow(k,i)=-10; %threw out because errors included neg.
                    end
                else
                    maglow(k,i)=-20; %didn't run because too small or descending
                end
            else
                maglow(k,i)=-30; %threw out because too few pts
            end
         end
    end
    
    
    tmp=nan(1,nx_geo);
    synth=zeros(size(dn));
    for i=1:length(dnr)
        %afid=find(dn>dnr(i));
        tmp(goodid)=mags(i,:);
        fwrite(fidmag(i),tmp,'real*4');
        tmp(goodid)=times(i,:);
        fwrite(fidt(i),tmp,'real*4');
        tmp(goodid)=maglow(i,:);
        fwrite(fidmagl(i),tmp,'real*4');
        tmp(goodid)=maghig(i,:);
        fwrite(fidmagh(i),tmp,'real*4');
        tmp(goodid)=timerr(i,:);
        fwrite(fidte(i),tmp,'real*4');
        tmp(goodid)=count(i,:);
        fwrite(fidn(i),tmp,'real*4');

    end
    
    res=-log(dat);
    resn=mean(res.^2,1,'omitnan');
    
    fwrite(fidres,resn,'real*4');
end
fclose('all');
