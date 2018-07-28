% unzips and renames data downloaded from earthdata ast14dem

% Paula Burgi
% written 06/13/18

clear 
sz = dir('03*'); 

% unzip 
for i = 1:length(sz)
    clear sf

sn = sz(i).name; 

if strfind(sn, '.zip')
   system(sprintf('unzip -o %s', sn)); 
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
met = cell2mat(y(4)); 

[~, dstr1] = system(sprintf('more %s | grep CALENDAR -A 2', met));
qidx  = strfind(dstr1,'"');
dstr  = dstr1(qidx(1)+1:qidx(2)-1); 

cd ..
system(sprintf('mv %s %s', sf, dstr)); 
cd(dstr); 
cd ..

end
