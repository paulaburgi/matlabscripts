% fix_lidar_dem_edges

close all

cdir = '/data/pmb229/other/OregonLidarDEM/temp/'; 
cd(cdir); 

x  = importdata('all_stitched_masked.dem.vrt');
l1 = x{1}; 
qf = strfind(l1, '"'); 
nx = str2num(l1(qf(1)+1:qf(2)-1)); 
ny = str2num(l1(qf(3)+1:qf(4)-1)); 


filename    = 'all_stitched_masked.dem'; 
fid         = fopen(filename,'r','native');
dem         = fread(fid,[nx,ny],'int16');
dem         = dem';  % dem(1:1000, 1:1000) is top left corner
% figure; imagesc(dem);
% caxis([0 2000]);

nnx = 100:1300; 
nny = 100:1300; 
dem2 = dem(nnx, nny); 
% figure; pcolor(dem2);
dem2 = dem; 

dem2v = dem2(:); 
sz    = size(dem2); 
dem3  = dem2; 
x     = round(linspace(1,length(dem2v), 100)');

g = zeros(size(dem2, 1), size(dem2, 2));
for i = 1:length(dem2v)
    [r,c] = ind2sub(sz, i); 
    if (r~=1) && (c~=1) && r~=sz(1) && c~=sz(2)
        neigh = [r+[-1;0;1;-1;1;-1;0;1] c+[-1;-1;-1;0;0;1;1;1]]; 
        ind   = sub2ind(sz,neigh(:,1),neigh(:,2));
        tmp   = dem2v(ind); %1:8
        tmp2  = dem2(r,c); %9
        nn    = length(find(tmp == 32767)); 
    if tmp2 == 32767
        if nn < 7 
%             g(r,c)    = 1; 
            nidx      = find(tmp ~= 32767); 
            %dem3(r,c) = mean(tmp(nidx));  
            dem3(r,c) = median(tmp(nidx));  
        end
    end
    if find(i == x)
        disp(find(i == x))
    end
    end
end

figure; 
subplot(2, 1, 1); 
imagesc(dem2); view(0,90); 
caxis([0 200]); 
subplot(2, 1, 2); 
imagesc(dem3); view(0,90); 
caxis([0 200]); 
linkaxes;


doprint = 0; 
if doprint == 1
    fid = fopen('OregonLidar_median.dem', 'w');
    fwrite(fid, dem3', 'int16');
    fclose(fid);
end











