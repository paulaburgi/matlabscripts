slcfol = '/data/pmb229/isce/p222f870/stack/merged/SLC/';
wfol   = '/data/pmb229/isce/p222f870/stack/looktest/';

cd(slcfol); 

% get names of all slcdirs
slcdir = dir('20*'); 
slcdir = {slcdir.name}'; 
nslc   = length(slcdir); 


for i=1 %:nslc-2
    sdir1 = cell2mat(slcdir(i));
    sdir2 = cell2mat(slcdir(i+1));
    sdir3 = cell2mat(slcdir(i+2));
    
    
    
end

























% slc 1
im = sqrt(-1); 
filename2    = [sdir1 '/' sdir1 '.slc']; %OR .FLAT.GEO 
fid          = fopen(filename2,'r','native');
[rmg2,count] = fread(fid,[nx*2,ny],'real*4');
status       = fclose(fid);
real1        = flipud((rmg2(1:2:nx*2,1:ny))');
imag1        = flipud((rmg2(2:2:nx*2,1:ny))');
mag          = abs(real1+im*imag1);
phs          = (angle(real1+im*imag1));

% slc 2
sdir1 = cell2mat(slcdir(2));
filename2    = [sdir1 '/' sdir1 '.slc']; %OR .FLAT.GEO 
fid          = fopen(filename2,'r','native');
[rmg2,count] = fread(fid,[nx*2,ny],'real*4');
status       = fclose(fid);
real2        = flipud((rmg2(1:2:nx*2,1:ny))');
imag2        = flipud((rmg2(2:2:nx*2,1:ny))');
mag          = abs(real2+im*imag2);
phs2         = (angle(real2+im*imag2));

% slc 3
sdir1 = cell2mat(slcdir(3));
filename2    = [sdir1 '/' sdir1 '.slc']; %OR .FLAT.GEO 
fid          = fopen(filename2,'r','native');
[rmg2,count] = fread(fid,[nx*2,ny],'real*4');
status       = fclose(fid);
real3        = flipud((rmg2(1:2:nx*2,1:ny))');
imag3        = flipud((rmg2(2:2:nx*2,1:ny))');
mag          = abs(real3+im*imag3);
phs3         = (angle(real3+im*imag3));
% 
% % interferograms
% [int1] =  circ_dist(phs,phs2);
% [int2] =  circ_dist(phs2,phs3);
% [int3] =  circ_dist(phs,phs3);
% 
% % phase triplets
% p1    = circ_dist(int1,-int2);
% ptrip = circ_dist(p1,int3);
% rl    = mag.*sin(p1); 
% ig    = mag.*cos(p1); 
% trp   = reshape([rl;ig], size(rl,1), []);


% get vrt template
tempvrt = [slcfol '20070713/20070713.slc.vrt'];
if ~exist(tempvrt, 'file')
    system(['cp ' tempvrt ' ' wfol]);
end
% write out interferograms
cd(wfol); 
for i = 1:3
    ri   = eval(['int' num2str(i)]);
    % make int
    rl    = mag.*sin(ri); 
    ig    = mag.*cos(ri); 
    ri   = reshape([ig;rl]', size(rl,1), []);
    % write int
    inti = ['int' num2str(i) '.int'];
    fid  = fopen([wfol inti], 'w'); 
    fwrite(fid, ri, 'real*4'); 
    fclose(fid);
    
    % make meta data
    system(['cp ' tempvrt ' ' inti '.vrt']); 
    system(['sed -i "s|20070713.slc|' inti '|g" ' inti '.vrt']); 
    system(['gdalbuildvrt ' inti '_2.vrt ' inti '.vrt']); 
    gdalbuildvrt2iscexml([inti '.vrt']); 
    system(['looks.py ' ]);
end

gdalbuildvrt2iscexml([wfol 'test2.int.vrt']);

xd = 1:3:nx; 
yd = 1:6:ny; %long dim 
figure; hold on; 
subplot(1,3,1); 
pcolor(phs(xd,yd)); shading flat; 
caxis([-3.14 3.14])
subplot(1,3,2); 
pcolor(phs2(xd,yd)); shading flat; 
caxis([-3.14 3.14])
subplot(1,3,3); 
pcolor(r1(yd,xd)'); shading flat; colormap jet; 
caxis([-3.14 3.14])













