% make_cctest_nonlinG.m

% make symbolic variables
z    = sym('z'); 
v    = sym('v'); 
tdem = sym('tdem'); 
svar  = [v; tdem; z];
% create function
G_m_ = double(G)*v + (tanh(mdc-tdem)'+1).*blgi*z.*0.5;
f    = G_m_ - double(intsiz); 
% create Jacobian 
J    = [diff(f, v) diff(f, tdem) diff(f, z)]; 
% initial guess & tolerance
var0  = [0 1 0];  
ep    = 1e-15; 

% LM method
[varf, k, Cm, X2, tdemall, zall, vall]   = LMLSQ_project(f, var0, J, ep, svar); 
varf