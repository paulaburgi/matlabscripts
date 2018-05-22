% run_ints_isce_allPols.m 

clear 

% parameters 
isceapp = 'insarApp.py'; % Alos=insarApp.py, sentinel=topsApp.py 

start = '--start=startup'; 
% start = ''; 

endend = '--end=filter'; 

% dostep = '--dostep=filter'; 

intdirs = dir; 
intdirs = {intdirs.name}; 

nints = []; 

for i=3:length(intdirs)
    % find number of polarization folders in int dir
    intdir = cell2mat(intdirs(i)); 
    cd(intdir)
    poldirs = dir; 
    poldirs = {poldirs.name}; 
    
    % cycle through polarization combos
    for j = 3:length(poldirs)
        poldir = cell2mat(poldirs(j)); 
        cd(poldir)
        disp(sprintf('\n\nrunning int %s, %s\n\n', intdir, poldir));
        system(sprintf('%s %s.xml --steps %s %s', isceapp, intdir, start, endend))
        cd ..
        nints = [nints; 1]; 
    end
    cd ..
    disp(sprintf('\n\n\n\n\n\n\n\n'))
end

sum(nints)
        