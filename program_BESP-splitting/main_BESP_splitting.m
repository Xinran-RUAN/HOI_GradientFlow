% main procedure: BESP-splitting
clear; para.dt=0.01; format long; % as dt -> 0, converge to the ground state
opt.isfixN=0; % ==0: stop until converge
              % ==1: stop after a fixed number of steps
% parameter setup (default)
if ~isfield(para,'NN'), para.NN=4; end;
if ~isfield(para,'R'), para.R=16; end;
if ~isfield(para,'beta'), para.beta=10; end;
if ~isfield(para,'delta'), para.delta=10; end;

para.dh=1/2^para.NN; para.M=2*para.R/para.dh; 
para.x=linspace(-para.R,para.R-para.dh,para.M)';
para.V=para.x.^2/2;
 
if ~isfield(opt,'Nmax'), opt.Nmax=30; end;

% local v.r.
x=para.x; V=para.V; R=para.R; NN=para.NN; 
M=para.M; dh=para.dh; dt=para.dt; 
beta=para.beta; delta=para.delta;
xi=(-M/2:M/2-1)'*pi/R;

% initial value setup
phi=1/pi^(1/4)*exp(-0.5*x.*x); 
rho=phi.^2;
fphi=myfft(phi); 
frho=myfft(rho);
E_i=R*sum(abs(fphi).^2.*xi.^2)+delta*R*sum(abs(frho).^2.*xi.^2)+dh*(sum(V.*phi.^2)+0.5*beta*sum(phi.^4));

phi_o=phi;

% get ground state
ts=tic;
[phi,E_list]=ground_solver_SP(phi,para,opt);
tsolve = toc(ts);
fprintf(' cpu: %.4f\n', tsolve);

E_list=[E_i,E_list];
