clear 
close all
im       = sqrt(-1);
intf='int_170308_170224'; 


nx=24660;
ny=2855; 

fid         = fopen(sprintf('%s/merged/filt_topophase.flat',intf),'r','native');

[rmg,count] = fread(fid,[nx*2,ny],'real*4');
  status      = fclose(fid);
  real        = flipud((rmg(1:2:nx*2,1:ny))');
  imag        = flipud((rmg(2:2:nx*2,1:ny))');
  mag         = abs(real+im*imag);
  phs         = (angle(real+im*imag));
  data        = (phs);
  
  
  %then mask 
  bb=[1600 2600 800 2200]; 
  data(1:1600,:)=0; 
  data(2600:end,:)=0; 
  data(:,1:800)=0; 
  data(:,2200:end)=0; 
  mag(data==0)=0;
  
  
  rl=(mag.*cos(data));
  img=(mag.*sin(data)); 
  
  real_diff=real-rl;
  figure
    pcolor(mag(1:10:end,1:100:end))
    shading flat
    
    rl=flipud(rl);
    img=flipud(img);
  
  %surf(data(1:5:end,1:5:end)); shading flat; colormap jet; set(gca, 'ydir', 'reverse'); %axis equal; 

  ap=zeros(size(data,1), size(data,2).*2); 
  ap(:,1:2:end)=rl; 
  ap(:,2:2:end)=img; 
  sAs=single(ap); 
  
  
  %then save 
  
% fid = fopen(sprintf('%s/merged/filt_topophase_masked.flat',intf), 'wb');
% fwrite(fid,sAs,'real*4');
% fclose(fid);

sAs = single([flipud(mag) flipud(data)]');
fid = fopen(sprintf('%s/merged/filt_topophase_masked.flat',intf),'wb');
fwrite(fid,sAs,'real*4');
fclose(fid);
  











  