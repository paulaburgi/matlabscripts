% loadhansen.m
function [t, lon, lat] = loadHansen(ly)

% ly = 'Hansen_GFC-2016-v1.4_lossyear_10N_100E.tif'; 

gt = geotiffread(ly);
t     = flipud(gt); 
[l,w] = size(t);

info = geotiffinfo(ly);
d    = info.PixelScale(1);
lon1=info.RefMatrix(3,1);
lat1=info.RefMatrix(3,2);


% define latitude and longitude start/stop
if strcmp(ly(end-4), 'E'); 
    lon = lon1 : d : lon1+(d*w)-d';
else
    lon = lon1-(d*w)+d : d : lon1';
end

if strcmp(ly(end-9), 'S'); 
    lat = lat1 : d : lat1+(d*w)-d';
else
    lat = lat1-(d*w)+d : d : lat1';
end



% % Plot
% londx      = 100;       latdx  = 100; 
% lonstr     = 1;         latstr = 1; 
% lonfin     = length(t); latfin = length(t); 
% %     latstr = lonstr; 
% %     latdx  = londx; 
% %     latfin = lonfin; 
%     
% lon2 = lon(lonstr:londx:lonfin); 
% lat2 = lat(latstr:latdx:latfin); 
% t2   = t(latstr:latdx:latfin, lonstr:londx:lonfin);
% 
% pcolor(lon2, lat2, t2); shading flat; 


















