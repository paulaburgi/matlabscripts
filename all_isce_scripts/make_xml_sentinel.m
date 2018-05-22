%Create ISCE topsApp files for a given set of dates, 
%needs a text file with the zip files date ordered 

%currently the script does not recognizes SDV images as master, so files
%must be sorted manually

clear all
clc
warning off

DEM='stitched_i2.dem';
unwrap='grass'; %unwrap='snaphu_mcf'; %icu/snaphu_mcf
filter='0.4';
bounding_box='37.4, 37.6, -121.81, -121.56';
region_interest='37.44,37.45,-121.76,-121.75';
posting='30';
swath='3';
orbdir='/data/ISCE/s1a/precise';
auxdir='/data/ISCE/s1a/aux_cal';
system('rm run_ints.sh');

system('ls -d1 S1*.zip > SAFE.txt')
cskh5=textread('SAFE.txt','%s');
%cskh5=textread('SAFE_t164.txt','%s');
cskh5=cellstr(cskh5);

%read SAR dates
[ncy,ncx]=size((cskh5));
for i=1:ncy
    aa=cskh5(i,:);
    aa=aa{1};
    yy(i)=str2num(aa(18:21));
    mm(i)=str2num(aa(22:23));
    dd(i)=str2num(aa(24:25));
end
tt=datenum(yy,mm,dd)';

%combine dates
a=combnk(tt,2); 
[masterSAR,index]=sort(a(:,2));
slaveSAR=a(index,1);

ind=masterSAR-slaveSAR;
masterSAR(ind>49) = []; %remove Bt> 48 days
slaveSAR(ind>49) = []; %remove Bt> 48 days
%masterSAR(ind>400) = []; %remove Bt> 1 year
%slaveSAR(ind>400) = []; %remove Bt> 1 year
num_ints=length(slaveSAR);


%sort SAFE file names, otherwise SDV -SSV modes will mess up the files
for i=1:length(cskh5)
    hh=cskh5(i);
    hh=char(hh);
    h(i,:)=hh(18:25);
end
dates=str2num(h(:,3:end));
h=str2num(h);
[h,h_pol_index]=sort(h);
dates=dates(h_pol_index);
cskh5=cskh5(h_pol_index); %reorder file names by date

%combine dates
a=combnk(dates,2);
[master,index]=sort(a(:,2));
slave=a(index,1);

%calculate time spans
m2=master+20000000;
s2=slave+20000000;
for i=1:length(m2)
    f2=num2str(m2(i));
    f1=num2str(s2(i));
    timespan(i) = datenum(str2num(f2(1:4)),str2num(f2(5:6)),str2num(f2(7:8))) - datenum(str2num(f1(1:4)),str2num(f1(5:6)),str2num(f1(7:8)));
    infos(i,:)= [str2num(f2) str2num(f1) timespan(i)];
end
% time_cutoff=0; %minimum time between acquisitions, change if an eruption is ocurring
% selected_dates=find(abs(timespan) >= time_cutoff);

%combine h5 names
ah5=combnk(cskh5,2);
%[masterh5,indexh5]=sort(ah5(:,2));
%slaveh5=ah5(indexh5,1);
masterh5=ah5(:,2);
slaveh5 =ah5(:,1);
masterh5=masterh5(index);
slaveh5 =slaveh5(index);

str=cellstr(strcat('int_',num2str(master),'_',num2str(slave),'.xml'));
fol=cellstr(strcat('int_',num2str(master),'_',num2str(slave)));

system('rm run_baseline.sh');


fid3 = fopen('run_ints.sh','a');
fprintf(fid3,'#!/bin/sh');
fprintf(fid3,'\n');
fclose(fid3);

%for jj=1:length(intproc)
for jj=1:(num_ints)
    
    %kk=char(intproc(jj));
    %kk=kk(1:17);
    km1=datestr(masterSAR(jj),'yy-mm-dd');
    kmm1=km1(1:2);
    kmm2=km1(4:5);
    kmm3=km1(7:8);
    ks1=datestr(slaveSAR(jj),'yy-mm-dd');
    kss1=ks1(1:2);
    kss2=ks1(4:5);
    kss3=ks1(7:8);
    kk=strcat('int_',strcat(kmm1,kmm2,kmm3),'_',strcat(kss1,kss2,kss3));
    fid2 = fopen('run_baseline.sh','a');
    fid3 = fopen('run_ints.sh','a');

    for ii=1:length(str)
        kk2=char(str(ii));
        kk2=kk2(1:17);

        if strcmp(kk2,kk)==1 %if xml == proc, create dir and xml file

            intfilename=char(str(ii));
            folder=char(fol(ii));
            mkdir(folder)
            fid = fopen(intfilename,'w');

            fprintf(fid,'<topsApp>\n');
            fprintf(fid,'	<component name="topsinsar">\n');
            fprintf(fid,'\n');
            fprintf(fid,'       <property name="Sensor Name">SENTINEL1</property>\n');
%             fprintf(fid,'       <component name="Dem">\n');
%             fprintf(fid,strcat('        <catalog>',DEM,'</catalog>\n'));
% %             fprintf(fid,'       </component>\n');
%             fprintf(fid,'\n');
%             fprintf(fid,'       <property name="doppler method">useDEFAULT</property>\n');
%             fprintf(fid,'\n');
%             fprintf(fid,strcat('       <property name="posting">',posting,'</property>\n'));

            fprintf(fid,'\n');
            fprintf(fid,'	<component name="master">\n');
            fprintf(fid,strcat('       <property name="safe">../',char(masterh5(ii)),'</property>\n'));
            fprintf(fid,'       <property name="output directory">master</property>\n');
            fprintf(fid,strcat('       <property name="orbit directory">',orbdir,'</property>\n'));
            fprintf(fid,strcat('       <property name="auxiliary data directory">',auxdir,'</property>\n'));
            fprintf(fid,strcat('	   <property name="swath number">',swath,'</property>\n'));
            fprintf(fid,strcat('	<property name="region of interest">[',region_interest,']</property>\n'));
            fprintf(fid,'	</component>\n');

            fprintf(fid,'\n');
            fprintf(fid,'	<component name="slave">\n');
            fprintf(fid,strcat('       <property name="safe">../',char(slaveh5(ii)),'</property>\n'));
            fprintf(fid,'       <property name="output directory">slave</property>\n');
            fprintf(fid,strcat('       <property name="orbit directory">',orbdir,'</property>\n'));
            fprintf(fid,strcat('       <property name="auxiliary data directory">',auxdir,'</property>\n'));
            fprintf(fid,strcat('	   <property name="swath number">',swath,'</property>\n'));
            fprintf(fid,strcat('	<property name="region of interest">[',region_interest,']</property>\n'));
            fprintf(fid,'	</component>\n');

            fprintf(fid,'\n');
            fprintf(fid,'	<property name="azimuth looks">1</property>\n');
            fprintf(fid,'	<property name="range looks">1</property>\n');
            fprintf(fid,strcat('	<property name="filter strength">',filter,'</property>\n'));
            fprintf(fid,'	<property name="do unwrap">True</property>\n');
            fprintf(fid,strcat('	<property name="unwrapper name">',unwrap,'</property>\n'));
            fprintf(fid,strcat('	<property name="demfilename">',DEM,'</property>\n'));
            
            
            fprintf(fid,'	<property name="geocode list">["merged/test2.unw"]</property>\n');
            fprintf(fid,strcat('	<property name="geocode bounding box">[',bounding_box,']</property>\n'));

            fprintf(fid,'\n');
            fprintf(fid,'</component>\n');
            fprintf(fid,'</topsApp>\n');

            fclose(fid);

            movefile(intfilename,folder);
            %system(strcat(['cp demLat_S42_S39_Lon_W074_W071.dem.wgs84 ',folder]))
            %system(strcat(['cp demLat_S42_S39_Lon_W074_W071.dem.wgs84.xml ',folder]))
            
%folder=char(fol(ii));
remove=strcat(['rm ',folder,'/*slc*']);
remove2=strcat(['rm ',folder,'/*raw*']);
string0=strcat(['cd ',folder]);
string=strcat(['topsApp.py ',folder, '.xml --steps --end=filter']);
disp(string)
%string=strcat(['insarApp.py ',folder, '.xml --steps --end=''preprocess''']);
string1=strcat(['cd ..']);

%fprintf(fid,remove);
%fprintf(fid,'\n');
%fprintf(fid,remove2);
%fprintf(fid,'\n');
fprintf(fid2,string0);
fprintf(fid2,'\n');
fprintf(fid2,string);
fprintf(fid2,'\n');
fprintf(fid2,string1);
fprintf(fid2,'\n');
%fprintf(fid,remove);
%fprintf(fid,'\n');
%fprintf(fid,remove2);
%fprintf(fid,'\n');
fclose(fid2);

%print file to run
string=strcat(['./run_isce.sh ',num2str(master(ii)),' ',num2str(slave(ii))]);
fprintf(fid3,string);
fprintf(fid3,'\n');
%end
fclose(fid3);


        else
        end

    end
end


% for ii=1:length(str)
% folder=char(fol(ii));
% remove=strcat(['rm ',folder,'/*slc*']);
% remove2=strcat(['rm ',folder,'/*raw*']);
% string0=strcat(['cd ',folder]);
% string=strcat(['insarApp.py ',folder,'/',folder, '.xml --steps --end=filter']);
% %string=strcat(['insarApp.py ',folder, '.xml --steps --end=''preprocess''']);
% string1=strcat(['cd ..']);
% 
% %fprintf(fid,remove);
% %fprintf(fid,'\n');
% %fprintf(fid,remove2);
% %fprintf(fid,'\n');
% fprintf(fid2,string0);
% fprintf(fid2,'\n');
% fprintf(fid2,string);
% fprintf(fid2,'\n');
% fprintf(fid2,string1);
% fprintf(fid2,'\n');
% %fprintf(fid,remove);
% %fprintf(fid,'\n');
% %fprintf(fid,remove2);
% %fprintf(fid,'\n');
% end
% fclose(fid);

%system('rm -r int* list* junk* mkdir* process* run_baseline.sh* *.ps')

%system('rm int*proc junk* list* process*')

system('chmod 777 run_ints.sh');