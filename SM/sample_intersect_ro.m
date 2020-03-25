dirs={'T101','T28'};
resampdir='T28_T101_resamp';
if(~exist(resampdir))
mkdir(resampdir);
end
chdir(resampdir);

for i=1:length(dirs)
	files=dir(['../' dirs{i} '/geo_VV/rel*.cor.geo']);
	for j=1:length(files)
		if(isempty(regexp(files(j).name,'_4r')))
		disp(files(j).name);
if(~exist(files(j).name))
	command=['gdalwarp -of ISCE -ot Float32 -te 54.830 17.781 55.742 20.598 -tr 0.00111 0.00111 -r med ../' dirs{i} '/geo_VV/' files(j).name ' ' files(j).name];	
system(command);
command=['gdalwarp -of vrt -ot Float32 -te 54.830 17.781 55.742 20.598 -tr 0.00111 0.00111 -r med ../' dirs{i} '/geo_VV/' files(j).name ' ' files(j).name '.vrt'];
system(command);
end
end
end
end

for i=1:length(dirs)
	 file=['../' dirs{i} '/geo_VV/c0.cor.geo'];
	 newfile=[dirs{i} '_c0.cor.geo'];
	 if(~exist(newfile))
command=['gdalwarp -of vrt -ot Float32 -te 54.830 17.781 55.742 20.598 -tr 0.00111 0.00111 -r med ' file ' ' newfile '.vrt'];
system(command)
command=['gdalwarp -of ISCE -ot Float32 -te 54.830 17.781 55.742 20.598 -tr 0.00111 0.00111 -r med ' file ' ' newfile ];
system(command)
end
end




%gdalwarp -of ISCE -ot Float32 -te 54.830 17.781 55.742 20.598 -tr 0.00111 0.00111 -r med ../T28/geo_VV/rel_20180531.cor.geo test3.geo

