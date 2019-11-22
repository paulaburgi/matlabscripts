%template
%note - for individual "real" plots, should plot at resolution of data if
%possible, not regridded.  This is for comparison between data types only.

%for each file type, loading function should extract lon/lat/data files and dn for
%each (potentially with time of day info)

%this code will then interpolate onto the lat/lon/date grid of interest,
%with different options.

cd /data/rlohman/HyperA/

datatype='CYGNSS';
%can be (add here) CYGNSS, ...
%place info about what is being loaded into the loading function.
%set bad data to nan within loading function.
pflag=0; %don't plot from within loading function
switch datatype
    case 'CYGNSS'
        data=feval(['read_' datatype],pflag);
end

dn    = [data.dn];
nd    = length(dn);
dateg = datenum(2018,01,01):datenum(2018,12,31);

scale=1;
switch scale
    case 0
        %this is area shown in bigger location map
        minlon=35;
        maxlon=65;
        minlat=10;
        maxlat=35;        
    case 1
        %this is zoom in region containing SAR track 130
        minlon=52;
        maxlon=60;
        minlat=15;
        maxlat=22;
end

%load simple coastline for cropping (need better one?)
load coast
coastlon=long;
coastlat=lat;
clear long lat


cd /data/pmb229/other/SM/cygnss

%target lat lon grids
dlat          = 0.01;
dlon          = 0.01;
ext           = 1; %extension for interpolation
[long ,latg]  = meshgrid([minlon:dlon:maxlon],[minlat:dlat:maxlat]);
lonboxbig=[minlon-ext minlon-ext maxlon+ext maxlon+ext];
latboxbig=[minlat-ext maxlat+ext maxlat+ext minlat-ext];

distcutoff=0.15; %degrees from closest obs -> nan
%regrid, within date, with cutoff
datacube=nan(11,size(long,1),size(long,2)); %datacube=nan(length(dateg),size(long,1),size(long,2));
for i=1:length(dateg)
    id1=find(abs(dateg(i)-dn)<1);
    if(isempty(id1))
        disp(['no data on date ' datestr(dateg(i))]);
    elseif(length(id1)==1)
        disp(datestr(dateg(i)));
        lon   = double(data(id1).lon);
        lat   = double(data(id1).lat);
        dat   = data(id1).data;
        
        id2   = inpolygon(lon,lat,lonboxbig,latboxbig);
        id3   = isfinite(dat);
        id    = and(id2,id3);
        [k,d] = dsearchn([lon(id) lat(id)],[long(:) latg(:)]);
        
        newg=griddata(lon(id),lat(id),dat(id),long,latg,'nearest');
        newg(d>distcutoff)=NaN;
        datacube(i,:,:)=newg;
      else
        disp(['too many dates? ' length(id1)])
    end
end  
    
outfile_name=[datatype '_' num2str(scale) '_loaded_d015_r001.mat'];
save(outfile_name,'datacube','dateg','long','latg');
   
%add stuff here to save as geotiffs, write csv file needed for geebam
%uploader?


close all; 
pcolor(long, latg, reshape(datacube(11,:,:), [size(datacube,2) size(datacube,3)])); shading flat;
hold on;
plot(coastlon, coastlat, 'k');



%% geotiff

dt = '7-23-2018'; 
dd = datenum(dt)-737060;
vv = reshape(datacube(dd,:,:), [size(datacube,2) size(datacube,3)]); 

d99 = vv; 
d99(isnan(d99)) = -99;
lonlim = [min(long(:)) max(long(:))]; 
latlim = [min(latg(:)) max(latg(:))]; 
s = size(d99);
R = georefcells(latlim, lonlim, s);
geotiffwrite(['cygnss_' dt '_c15.tif'], d99, R);































