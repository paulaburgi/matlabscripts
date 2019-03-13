% run_ints in roipace
% Paula Burgi, 2017

clear 

addpath('/home/pmb229/matlab');

ds = 'ds_bl-1000m_dl-500m.mat'; 
DEMtype = 'NED'; 
startFrom = 'raw'; 
endAt = 'done'; 

ds_exist = exist(ds); 

sz = dir('int_*proc'); 

ds2 = struct(); 
if ds_exist == 2
    load(ds); 
    ds1 = ds.ds1; 
    ds2 = ds.ds2; 
else
    ds1=[]; ds2=[];
    for i=1:length(sz)
        fprintf(repmat('\n',1,20));
        sn = sz(i).name; 
        dd1 = sn(5:10);
        dd2 = sn(12:17); 
        ds1 = [ds1; dd1]; 
        ds2 = [ds2; dd2];
    end
end

for i=1:length(ds1)
    fprintf(repmat('\n',1,20));
    d1 = ds1(i,:);
    d2 = ds2(i,:);
    folname = sprintf('int_%s_%s_%s', d1, d2, DEMtype); 
    pname = sprintf('%s.proc', folname); 
    curd = pwd; 
    if exist(folname, 'dir'); 
        % start does not equal raw
        if strcmp(startFrom,'full_res') | strcmp(startFrom,'begin_filt')
%             cd(folname);
%             system('rm -rf *unw *HDR* *SIM* Geo/ radar_4rlks.hgt* geo*'); 
%             cd(curd); 
%             system(sprintf('process_2pass.pl %s %s %s', pname, startFrom, endAt)); 
            cd(folname);
            system(sprintf('touch %s-%s.cor', d1, d2));         
            cd(curd); 
            system(sprintf('process_2pass.pl %s %s %s', pname, startFrom, endAt)); 
        elseif strcmp(startFrom, 'offsets')
            cd(folname);
            system(sprintf('touch %s.slc', d1)); 
            system(sprintf('touch %s.slc', d2));
            system(sprintf('rm %s-%s.int', d1, d2)); 
            system(sprintf('rm %s-%s.amp', d1, d2));
            cd(curd); 
            system(sprintf('process_2pass.pl %s %s %s', pname, startFrom, endAt)); 
        elseif strcmp(startFrom, 'unwrapped')
            cd(folname);
            system('touch phase_var_*_*rlks.msk.rsc'); 
            cd(curd); 
            system(sprintf('process_2pass.pl %s %s %s', pname, startFrom, endAt)); 
        elseif strcmp(startFrom, 'raw') && exist(folname, 'dir')
            fprintf('Skipping %s, int folder exists\n', pname);  
        end
        
    elseif strcmp(startFrom, 'raw') && exist(folname, 'dir')
        fprintf('Skipping %s, int folder does not exist\n', pname); 
    elseif strcmp(startFrom, 'raw') && exist(folname, 'dir')
        system(sprintf('process_2pass.pl %s %s offsets', pname, startFrom));
            cd(folname)
            newfitoff(d1, d2);
            cd(curd)
            close all
            system(sprintf('process_2pass.pl %s %s %s', pname, 'offsets', endAt));
    end
end
% else
%     for i=1:length(sz)
%         fprintf(repmat('\n',1,20));
%         sn = sz(i).name; 
%         d1 = sn(5:10);
%         d2 = sn(12:17); 
%         folname = sprintf('int_%s_%s_%s', d1, d2, DEMtype); 
%         if strcmp(startFrom,'full_res') && exist(folname, 'dir')
%             curd = pwd; 
%             cd(folname);
%             system('rm -rf *unw *HDR* *SIM* Geo/ radar_4rlks.hgt* geo*'); 
%             cd(curd); 
%             system(sprintf('process_2pass.pl %s %s %s', sn, startFrom, endAt)); 
%         elseif strcmp(startFrom,'full_res') && ~exist(folname, 'dir')
%             fprintf('Skipping %s, int folder does not exist\n', sn); 
%         elseif strcmp(startFrom, 'fill these in')
%             %do something
%         elseif strcmp(startFrom, 'raw') && ~exist(folname, 'dir')
%             system(sprintf('process_2pass.pl %s %s offsets', sn, startFrom));
%             curd = pwd; 
%             cd(folname)
%             newfitoff(d1, d2);
%             cd(curd)
%             close all
%             system(sprintf('process_2pass.pl %s %s %s', sn, 'offsets', endAt));
%         elseif strcmp(startFrom, 'raw') && exist(folname, 'dir')
%             fprintf('Skipping %s, int folder exists\n', sn); 
%         end
%     end