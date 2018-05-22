function newfitoff(date1,date2)

name1=[date1 '-' date2 '_ampcor.off'];
name2=[date1 '-' date2 '_cull.off'];


a=load(name1);
b=load(name2);

figure
plot(a(:,1),a(:,4),'.')
hold on
plot(b(:,1),b(:,4),'r.')
legend('alloff','ROI_PAC cull');
xlabel('range direction');
ylabel('azimuth offset');

x=a(:,1);
y=a(:,3);
dx=a(:,2);
dy=a(:,4);
stds=a(:,6);

figure
subplot(1,4,1)
scatter(x,y,10,stds,'filled')
axis tight
colorbar('h')

subplot(1,4,2)
scatter(x,y,10,a(:,7),'filled')
axis tight
colorbar('h')

subplot(1,4,3)
scatter(x,y,10,dy,'filled')
axis tight
colorbar('h')

subplot(1,4,4)
scatter(x,y,10,a(:,5),'filled')
axis tight
colorbar('h')

firstgood=find(stds<.03);
if(length(firstgood)<.05*length(x))
disp('too few points retained in initial cutoff')
    firstgood=find(stds<.2);
end

np=length(x);
id=1:np;
id=firstgood;
mx=x(id);
my=y(id);
d=dy(id);
all=a(id,:);
for i=1:2
    G=[ones(size(mx)) mx my mx.*my mx.^2 my.^2];


    [mod,resn,res]=lsqlin(G,d);

    synth=G*mod;
    cutoff=std(res);
    id=find(abs(res)>cutoff);
    id2=find(abs(res)<=cutoff);

    
    figure
    subplot(1,4,1)
    scatter(x,y,10,dy,'filled')
    hold on
    plot(mx,my,'k.','markersize',1')
    %axis image
    title('all data and data used in this step')
     colorbar('h')
    axis tight   
    subplot(1,4,2)
    scatter(mx(id2),my(id2),10,synth(id2),'filled')
    %axis image
    title('synthetic')
     colorbar('h')
    axis tight   
    subplot(1,4,3)
    scatter(mx(id2),my(id2),10,d(id2),'filled')
    %axis image
    title('picked')
 colorbar('h')
axis tight
 subplot(1,4,4)
    scatter(mx(id2),my(id2),10,d(id2)-synth(id2),'filled')
    colorbar('h')
    axis tight
    mx=mx(id2);
    my=my(id2);
    d=d(id2);
    all=all(id2,:);
end

   Gall=[ones(size(x)) x y x.*y x.^2 y.^2];


    %[mod,resn,res]=lsqlin(G,d);

    synthall=Gall*mod;
    resall=dy-synthall;
    %id=find(abs(res)>cutoff);
    id3=find(abs(resall)<=cutoff);


figure
subplot(1,2,1)
plot(a(:,1),a(:,4),'.')
hold on
plot(mx,d,'r.')
plot(x(id3),dy(id3),'w.','markersize',1)
legend('alloff','culled');
xlabel('range direction');
ylabel('azimuth offset');

subplot(1,2,2)
plot(x,y,'.')
hold on
plot(mx,my,'ro')
plot(x(id3),y(id3),'g.')
legend('all','picked')
title('points picked in x,y')

system(['mv ' name2 ' oldcull.off']);

fid=fopen(name2,'w');

form='%8d%10.3f%8d%12.3f%11.5f%11.6f%11.6f%11.6f\n';
for i=1:length(id2)
    fprintf(fid,form,all(i,:));
end
fclose(fid);

