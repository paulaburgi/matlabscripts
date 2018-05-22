function int2kml_phase_only(infile, outfile)
% e.g. int2kml_phase_only('filt_topophase.flat.geo', 'sumatra_070703_070818_phs.png', 1513, 768)
% infile='geo_rates_2.unw';
% outfile='ENVI_T485_2889';
mode=2;
wraprate=6.2832;
scale=100;
%Usage: unw2png.pl imagefile [outfile] [mode] [scale] [justpng] 
%Function:  runs mdx based on imagefile.rsc, then generates ppm file with 
%   transparency based on pixel value in upper left corner
%mode: plot amp (1), phs (2), or both (3)  default = 2
%justpng: 1 means mdx has already been run and existing out.ppm is used
%scale in \%, default=100;
% WRITTEN BY KYLE MURRAY

% set_params
% load(ts_paramfile);
% if strcmp(sat,'S1A')
%     nx=ints(id).width;
%     ny=ints(id).length;
% else
%     [nx,ny,lambda] = load_rscs(infile,'WIDTH','FILE_LENGTH','WAVELENGTH');
% 
% end

lamda = 0.236; % for ALOS stripmap

% find nx, ny
x=importdata([infile '.vrt']);
l1 = x{1}; 
qf = strfind(l1, '"'); 
nx = str2num(l1(qf(1)+1:qf(2)-1)); 
ny = str2num(l1(qf(3)+1:qf(4)-1)); 


comman = ['cpx2mag_phs ' infile ' mag phs ' num2str(nx)];
system(comman);

   if mode==1
        comman=['mdx mag -r4 -ponly '  num2str(nx)];
        system(comman);
   else
       comman=['mdx phs -r4 -cmap cmy -wrap ' num2str(wraprate) '  -ponly ' num2str(nx)];
       system(comman);
   end
       
!rm mag phs

% read first pixel
comman = ['convert out.ppm -crop 1x1+0+0 txt:->tmp'];
system(comman);

color  = 1;
!rm tmp

comman = ['convert out.ppm -resize ' num2str(scale) '%' ' -transparent cyan ' outfile];

system(comman);
% print STDERR "got here!\n";
comman =  ['importIMGtoKML.pl ' outfile ' ' infile ' ' outfile '.kml'];
system(comman);


% copy infile and run normal mdx.py -kml, to get correct cordinates for
% .kml file created in the prior part of this script. 
infile2 = [infile(1:end-4) '2' infile(end-3:end)];
system(['cp ' infile ' ' infile2])
system(['cp ' infile '.xml ' infile2 '.xml'])
system(['mdx.py ' infile2 ' -kml mattemp.kml']);
mattemp = importdata('mattemp.kml'); 

if exist([outfile '.kml'], 'file')
    tmp = fileread([outfile '.kml']); tmp2 = tmp; 
        tmp2 = strrep(tmp2, '<north></north>', strtrim(mattemp{contains(mattemp, 'north')})); 
        tmp2 = strrep(tmp2, '<south>0</south>', strtrim(mattemp{contains(mattemp, 'south')})); 
        tmp2 = strrep(tmp2, '<east>0</east>', strtrim(mattemp{contains(mattemp, 'east')})); 
        tmp2 = strrep(tmp2, '<west></west>', strtrim(mattemp{contains(mattemp, 'west')})); 
        
    fid = fopen([outfile '.kml'], 'wt'); 
    fwrite(fid, tmp2); 
    fclose(fid); 


end

system('rm mattemp.kml'); 
system(['rm ' infile2 '.png']);
system(['rm ' infile2]); 
system(['rm ' infile2 '.xml']); 

