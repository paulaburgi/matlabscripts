% make raw SLCs for ALOS data for roipac processing
% unzips downloaded data, changes name to date, and runs make_raw_alos.pl

% Paula Burgi
% written 10/17/17

clear 
sz = dir('ALPS*'); 

% unzip 
for i = 1:length(sz)
    clear sf

sn = sz(i).name; 

if strfind(sn, '.zip')
   system(sprintf('unzip %s', sn)); 
   if exist('zipdata')
       system(sprintf('mv %s zipdata',sn))
   else
   system(sprintf('rm %s', sn));
   end
   sf = sn(1:end-4);
else
    sf = sn;
end

cd(sf)
x = dir; 
y = {x.name};
z = cell2mat(y); 

if strfind(z,'HH') & strfind(z,'HV')
    p = 'FBD2FBS'; 
    else
    p = 'NO'; 
end

[~, dstr1] = system(sprintf('more %s.l0.workreport | grep Img_SceneCenterDateTime',sf(1:end-5)));
dstr2 = dstr1(strfind(dstr1,'='):end); 
dstr = dstr2(6:11); 

cd ..

system(sprintf('mv %s %s', sf, dstr)); 

cd(dstr); 

system(sprintf('make_raw_alos.pl IMG %s %s', dstr, p)); 

cd ..

end
