% main procedure: BEFD-splitting
clc;clear; para.dt=1e-2;
opt.isfixN=0; % ==0: stop until converge
              % ==1: stop after a fixed number of steps
% parameter setup (default)
if ~isfield(para,'NN'), para.NN=3; end;
if ~isfield(para,'R'), para.R=16; end;
if ~isfield(para,'beta'), para.beta=10; end;
if ~isfield(para,'delta'), para.delta=10; end;
if ~isfield(para,'dt'), para.dt=1e-3; end;

para.dh=1/2^para.NN; para.M=2*para.R/para.dh; 
para.x=linspace(-para.R,para.R,para.M+1)';
para.V=para.x.^2/2;

if ~isfield(opt,'isfixN'), opt.isfixN=0; end;
if ~isfield(opt,'Nmax'), opt.Nmax=50; end;
if ~isfield(opt,'state'), opt.state=0; end;

% local v.r.
x=para.x; V=para.V; 
NN=para.NN; 
M=para.M;
dh=para.dh; 
beta=para.beta; delta=para.delta;


% initial value setup
if opt.state==0
    phi=1/pi^(1/4)*exp(-0.5*x.*x);
else
    phi=sqrt(2)*1/pi^(1/4)*x.*exp(-0.5*x.*x);
end
dphi=(phi(3:M+1)-phi(1:M-1))/(2*dh);
E_i=dh*(sum(0.5*dphi.^2)+sum(V.*phi.^2)+0.5*beta*sum(phi.^4)+sum(2*delta*phi(2:M).^2.*dphi.^2));


% get ground state
[phi,E_list]=ground_solver(phi,para,opt);

E_list=[E_i,E_list];

%% save fig, data, ...
fname = strcat('GPE-FD1d-Bet-',num2str(beta),'-Del-',num2str(delta),'-NN-',num2str(NN));
E=E_list(end); save(strcat(fname,'.mat'),'E','phi','x');
