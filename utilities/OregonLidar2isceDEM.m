% OregonLidar2isceDEM.m

cdir = pwd; 

% unzip any zipped lidar folders
zdir = dir('LDQ*zip'); 
zdir = {zdir.name};

if ~isempty(zdir) 
    for i = 1:length(zdir)
        zdiri = cell2mat(zdir(i)); 
        system(['unzip ' zdiri]); 
        system(['rm ' zdiri]);
    end
end

Ldir = dir('LDQ*');
Ldir = {Ldir.name};

Ll = length(Ldir); 

for i = 1:Ll 
    Ldiri = cell2mat(Ldir(i)); 
    indir = ['be' Ldiri(5:9) lower(Ldiri(10)) Ldiri(11)]; 
    
    if exist(['downlooked/' indir '.dem.xml'], 'file'); 
        % do nothing
    else
        
    cd(Ldiri)
    
    % check if there are more options as downloads finish... 
    if exist('2015_OLC_Upper Umpqua 3DEP', 'dir')
        adir = ['2015_OLC_Upper Umpqua 3DEP/Bare_Earth/' indir '/'];
    else
        adir = ['2015_OLC_Lane County/Bare_Earth/' indir '/'];
    end
    cd(adir);
    
    % !!!!!!!!!1add Work Log work flow here!!!!!!!!!!!!
    
    
    cd(cdir);
    end
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    