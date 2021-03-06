% already have Aster DEM tif's in 1 folder. 
% Make 1 master csv for all files. 

clear 
close all

%% define files, folders 
% cascadia
pf_fol  = '/data/pmb229/isce/p222f870/DEMs/Aster/'; 
geefol  = [pf_fol 'AsterDEM_forGEE/'];

cd(pf_fol); 
demdirs = dir('20*'); 
demdirs = {demdirs.name}; 
ix = 0; 
dstra = [];

fid = fopen([geefol 'meta_all.csv'], 'wt'); 
fprintf(fid, 'id_no,idx,date,year\n'); 

for i = 1:length(demdirs)
    cd(cell2mat(demdirs(i)))
    
    % get info from meta files
    x          = dir; 
    y          = {x.name}; 
    met        = cell2mat(y(contains(y, 'met')));
    [~, dstr1] = system(sprintf('more %s | grep CALENDAR -A 2', met));
    qidx       = strfind(dstr1,'"');
    dstr       = dstr1(qidx(1)+1:qidx(2)-1);
    year       = dstr1(qidx(1)+1:qidx(1)+4);
        dstra = [dstra; dstr];
    
    % write csv
    fprintf(fid, [met(1:end-8) ',']); 
    fprintf(fid, [num2str(ix) ',']); 
    fprintf(fid, [dstr ',']); 
    fprintf(fid, [year '\n']); 
    
    ix = ix+1; 
    cd ..
end   
    
fclose(fid); 









