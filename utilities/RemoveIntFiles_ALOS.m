% Remove unecessary files after processing with InsarApp
% RemoveIntFiles_ALOS.txt is in /home/pmb229/matlab/utilities/

% d=importdata('RemoveIntFiles_ALOS.txt'); 
d=importdata('RemoveIntFiles_ALOS_extra2.txt'); 

% directory with int folders
% dr = '/data/pmb229/isce/p222f870/mostcombos/'; 
dr = '/data/pmb229/isce/p222f870/HVcombos/'; 
cd(dr); 
intdir = dir('int_*'); 

for j = 1:length(intdir)
    cd(intdir(j).name)
    
    for i = 1:length(d)
        di = cell2mat(d(i)); 
        system(['rm ' di ' ' di '.xml ' di '.vrt']); 
    end
    
%     dts = intdir(j).name; 
%     system(['rm ' dts(5:10) ' ' dts(12:end) ' 20' dts(5:10) '.slc.geo*' '  20'  dts(12:end) '.slc.geo*']); 
    
    cd ..
end
