pffol       = '/data/pmb229/isce/p222f870/'; 
intfol      = [pffol 'mostcombos/']; 
datafol     = [pffol 'data/']; 

cd(intfol); 

datdirs  = dir([datafol '20*']); 
datdirs2 = cell2mat({datdirs.name}'); datdirs = [];
for i=1:length(datdirs2)
    datdirs = [datdirs; str2double(datdirs2(i,:))];
end
intdirs  = dir('int_*'); 
intdirs  = {intdirs.name}';

intdates = []; 
for i = 1:length(intdirs); 
    idir = cell2mat(intdirs(i,:)); 
    d1   = str2double(['20' idir(5:10)]); 
    d2   = str2double(['20' idir(12:end)]); 
    intdates = [intdates; d1 d2]; 
end


for i = 1:length(datdirs)
    didx = find(datdirs(i) == intdates(:,1));
    iss = 'master';
    if isempty(didx)
        didx = find(intdates(end) == intdates(:,2));
        iss = 'slave'; 
        idr = cell2mat(intdirs(didx(end))); 
    else
        idr = cell2mat(intdirs(didx(1))); 
    end
    cd(idr); 
    
    dldir = dir([iss '*dl*']);
    if isempty(dldir)
        system(['looks.py -i ' iss '.slc -o ' iss '_dl.slc -r 3 -a 6']); 
    end
    
    gcdir = dir([iss '_dl*geo']);
    if isempty(gcdir)
        system(['insarApp.py ' idr '.xml --steps --dostep=geocode']); 
    end
    
    gcslc = [iss '_dl.slc.geo']; 
    if strcmp(iss, 'master'); 
        nf = ['20' idr(5:10) '.slc.geo']; 
        system(['cp ' gcslc ' ' nf]);  
        system(['cp ' gcslc '.vrt ' nf '.vrt']);  
        system(['cp ' gcslc '.xml ' nf '.xml']);  
    else
        nf = ['20' idr(12:end) '.slc.geo']; 
        system(['cp ' gcslc ' ' nf]);  
        system(['cp ' gcslc '.vrt ' nf '.vrt']);  
        system(['cp ' gcslc '.xml ' nf '.xml']);  
    end
    
    system(['sed -i s/' gcslc '/' nf '/g ' nf '.vrt']);
    system(['sed -i s/' gcslc(1:end-4) '/' nf(1:end-4) '/g ' nf '.xml']);
    system('cp  20*.slc.geo* ../geocodedSLCs/'); 
    
    
    cd ..
end

    






