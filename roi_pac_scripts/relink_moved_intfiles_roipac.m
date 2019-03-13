% thisvSIM_4rlks.hgt script re-links all the files within a roipac interferogram
% directory whose links are broken due to moving around files 

clear

% in Geo directory
f1 = 'geomap_4rlks.trans'; 
f2 = 'geomap_4rlks.trans.rsc'; 
f3 = 'geomap_4rlks.trans.rsc.hst'; 
% in SIM directory 
f4 = 'SIM_16rlks.hgt'; 
f5 = 'SIM_16rlks.hgt.rsc'; 
f6 = 'SIM_16rlks.hgt.rsc.hst'; 
f7 = 'SIM_4rlks.hgt'; 
f8 = 'SIM_4rlks.hgt.rsc'; 
f9 = 'SIM_4rlks.hgt.rsc.hst'; 
f10 = 'Geo/geomap.orrm'; 

% extract dates out of path string 
% slashidx = find(cdir =='/'); 
% 
% % if folder is named 'int_d1_d2_whatever'
% fld = cdir(slashidx(end)+5:end); 
% usidx = find(fld =='_'); 
% fld(usidx) = ' '; 
% dates = sscanf(fld, '%g', 2); 


alldir=dir; 
namedir={alldir.name}; 
subdirs = namedir([alldir.isdir]);
ud = pwd; 

for k=1:length(subdirs); 
    d = subdirs{k}; 
    if exist(sprintf('%s/log', d))
        cd(d)
        cdir = pwd;
       for i=1:3
        f = eval(sprintf('f%g', i));
        nf = sprintf('%s/Geo/%s', cdir, f); 
        system(sprintf('ln -sf %s %s', nf, f)); 
    end
    for i=4:10
        f = eval(sprintf('f%g', i)); 
        if i == 10
            nf = sprintf('%s/SIM/%s', cdir, 'SIM.orrm'); 
        else
        nf = sprintf('%s/SIM/%s', cdir, f); 
        end
        system(sprintf('ln -sf %s %s', nf, f)); 
    end
    else 
        %pass 
    end
    cd(ud)
end
    



    
    
    
    
    
    
    
    
    
    
    
    
    