%startup.m

%addpath( '/home/pmb229/matlab/ALLPS/resamp'); 
%addpath /home/pmb229/matlab/ALLPS/
addpath /home/rlohman/insar/MY_SCR/
addpath /home/rlohman/insar/INT_SCR/

setenv('GFORTRAN_STDIN_UNIT','5');
setenv('GFORTRAN_STDOUT_UNIT','6');
setenv('GFORTRAN_STDERR_UNIT','0');
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin:/opt/local/lib:');

setenv('LD_LIBRARY_PATH','/usr/lib:/usr/local/lib:/usr/local/lib/x86_64-linux-gnu');