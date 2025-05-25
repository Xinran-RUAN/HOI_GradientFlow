% main procedure: BEFD-splitting-multiscale
clc; clear; tic;
para.dt=1e-2; opt.NNmax=7;
opt.isfixN=0; % ==0: stop until converge
              % ==1: stop after a fixed number of steps
% parameter setup (default)
if ~isfield(para,'NN'), para.NN=1; end
if ~isfield(para,'R'), para.R=16; end
if ~isfield(para,'beta'), para.beta=10; end
if ~isfield(para,'delta'), para.delta=1e4; end

NN=para.NN;

while(para.NN<opt.NNmax+1)
    para.dh=1/2^para.NN; para.M=2*para.R/para.dh; 
    para.x=linspace(-para.R,para.R,para.M+1)';
    para.V=para.x.^2/2;
    % local v.r.
    x=para.x; V=para.V; M=para.M; dh=para.dh;  
    beta=para.beta; delta=para.delta;
    if (para.NN==NN)
        % initial value setup
        phi=1/pi^(1/4)*exp(-0.5*x.*x);
    else
        phi=finer(phi,para);
    end
    % get ground state
    [phi,~]=ground_solver(phi,para,opt);
    dphi=(phi(3:M+1)-phi(1:M-1))/(2*dh);
    E=dh*(sum(0.5*dphi.^2)+sum(V.*phi.^2)+0.5*beta*sum(phi.^4)+sum(2*delta*phi(2:M).^2.*dphi.^2));
    para.NN=para.NN+1;    
end

toc

plot(para.x,phi)
disp(E)
% % save .mat
% fname = strcat('GPE-FD1d-Bet-',num2str(beta),'-Del-',num2str(delta),'-NN-',num2str(para.NN-1));
% save(strcat(fname,'.mat'),'E','phi','x');



