params
dem='/data/pmb229/isce/p222f870/Sentinel_p13f445/dem/all_stitched_masked.dem';
bbox='43.8 44.1 -123.5 -123.2';
%T130
dem='/data/rlohman/Sentinel/Saudi/dem/demLat_N17_N22_Lon_E054_E058.dem.wgs84';
bbox='17.75 21.93 54.4 57.91';
geodir=['geo_' pol];

if(~exist(geodir,'dir'))
    mkdir(geodir)
end

resdir=['results_dates_' pol '/'];

filelist=dir([resdir '*cor']);
nfiles=length(filelist)
for i=1:nfiles
    file=filelist(i).name;
    if(or(exist([geodir '/' file '.geo'],'file'),exist([resdir file '.geo'],'file')))
        disp([file ' already geocoded'])
    else
        geocoder(file,pol,resdir,bbox,dem,geodir);
    end
end


if(exist([geodir '/rows.geo'],'file'))
    disp('rows already geocoded')
else
    fid=fopen([resdir 'rows'],'w');
    for i=1:newny
        fwrite(fid,i*ones(1,newnx),'real*4');
    end
    fclose(fid);
    geocoder('rows','',resdir,bbox,dem,geodir);
    movefile([resdir 'rows.geo'],[geodir '/rows.geo']);
end


if(exist([geodir '/cols.geo'],'file'))
    disp('rows already geocoded')
else
    fid=fopen([resdir 'cols'],'w');
    for i=1:newny
        fwrite(fid,1:newnx,'real*4');
    end
    fclose(fid);
    geocoder('cols','',resdir,bbox,dem,geodir);
    movefile([resdir 'cols.geo'],[geodir '/cols.geo']);
end





    