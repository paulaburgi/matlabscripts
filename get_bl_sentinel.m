d1dir = dir('int_170224*'); 
d1dir = {d1dir.name}; 

bl_all = 0;
dn_all = datenum(['170224'], 'yymmdd');
for i = 1:length(d1dir)
    di = cell2mat(d1dir(i)); 
    dn_all = [dn_all; datenum(di(12:end), 'yymmdd')];
    
    cd(di); 
    [~,t]  = system('grep Bperp isce.log'); 
    eidx = strfind(t, '=');
    bl = str2double(t(eidx(1)+1:eidx(1)+9));
    bl_all = [bl_all; bl];
    cd ..
end

dcom = combnk(dn_all,2);
bcom = combnk(bl_all,2);

bl_lim = 1000; 
dt_lim = 25;

diffd = abs(diff(dcom'));
bc = find(diffd<dt_lim);
dcom_tb = dcom(bc,:); 
bcom_tb = bcom(bc,1)-bcom(bc,2); 

bl=struct; 
bl.dn_all = dn_all; 
bl.bl_all = bl_all; 
dc=struct; 
dc.dateCombos = dcom; 
dc.bl  = bcom(:,2) - bcom(:,1); 
dc.good_cor_idx = bc'; 

save('cctest_sentinel_bl', 'bl'); 
save('cctest_sentinel_dc', 'dc'); 

keyboard