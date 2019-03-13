% replace_string_xml.m

d = dir('int_*'); 
d = {d.name};

old = '43.77 43.91 -123.44 -123.30';
new = '43.58 44.38 -123.76 -122.96';

old = '"filt_topophase.flat", "los.rdr", "topophase.cor", "topophase.flat", "phsig.cor"'; 
new = '"filt_topophase.flat", "topophase.flat"';


for i = 1:length(d)
    di = cell2mat(d(i));
    
    % text with spaces
    %    system(['sed -i "s:' old ':' new ':g" ' di '/' di '.xml']);
    
    % text with spaces and quotes
        system(['sed -i ''s:' old ':' new ':g'' ' di '/' di '.xml']);
        
end
