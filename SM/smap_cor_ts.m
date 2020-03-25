cd /data/rlohman/Sentinel/Saudi/T130_T28_resamp/

% directory
d = dir('rel_*.cor.geo'); 
d = {d.name}; 

% get dates
for i = 1:length(d)
	n1    = cell2mat(d(i)); 
	dn(i) = datenum(n1(5:12), 'yyyymmdd'); 
end
% storm date
sd = find(dn == datenum('20180526', 'yyyymmdd')); 

% size
nx = 822;
ny = 2538;

% pixels of interest
xpt = [263 309];
ypt = [1128 1092];

% pull out rel cor values at pixels of interest
for j = 1:length(xpt);
for i = 1:length(d)
	ni = cell2mat(d(i)); 
	fid = fopen(ni, 'r');
	fseek(fid, (nx*(ypt(j))+xpt(j))*4, -1);
	p = fread(fid, 1, 'real*4'); 
	fclose(fid); 
	v(i,j) = p;
end
end

% get timespan and mag of fit
fid  = fopen('/data/rlohman/Sentinel/Saudi/results_TS_T130_T28_VV/20180524.time0');
fid2 = fopen('/data/rlohman/Sentinel/Saudi/results_TS_T130_T28_VV/20180524.mag0');
t0   = fread(fid, [nx ny], 'real*4')'; fclose(fid);
m0   = fread(fid2, [nx ny], 'real*4')'; fclose(fid2);
a    = -log(diag(m0(ypt, xpt))); 
b    = 1./diag(t0(ypt, xpt));
x    = [1: 0.1:100]; 
for i = 1:length(a)
	f(i,:) = a(i)*exp(-b(i)*x);
end
f1    = a(1)*exp(-b(1)*x);
f2    = a(2)*exp(-b(2)*x); 


close all; 
 figure; hold on; box on;  
 cm = [0.8 0.5 0.1; 0.8 0.1 0.1]; 
 ms = 15; 
 h1 = plot(x, 1-f(1,:), '-', 'linewidth', 2, 'color', cm(1,:)); 
 h2 = plot(x, 1-f(2,:), '-', 'linewidth', 2, 'color', cm(2,:));
 plot(dn(1:sd-1)-dn(sd)+2, v(1:sd-1,1), '.', 'markersize', ms, 'linewidth', 1, 'color', cm(1,:));
 plot(dn(1:sd-1)-dn(sd)+2, v(1:sd-1,2), '.', 'markersize', ms, 'linewidth', 1, 'color', cm(2,:));
 plot(dn(sd:end)-dn(sd)+2, v(sd:end,1), '.', 'markersize', ms, 'linewidth', 1, 'color', cm(1,:));
 plot(dn(sd:end)-dn(sd)+2, v(sd:end,2), '.', 'markersize', ms, 'linewidth', 1, 'color', cm(2,:));
 xlim([-30 90]); 
 %xlim([datenum('20180501', 'yyyymmdd') datenum('20180815', 'yyyymmdd')]);
 %datetick('keeplimits'); 
ylim([0 1.05]); 
 grid minor; 
xlabel('Days since storm'); 
ylabel('relative coherence'); 
lgd=legend([h1, h2], ['pixel ' num2str(xpt(1)) ', ' num2str(ypt(1))], ['pixel ' num2str(xpt(2)) ', ' num2str(ypt(2))], 'location', 'southeast');
  
