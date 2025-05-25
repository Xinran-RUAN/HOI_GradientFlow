% main procedure: BEFD-splitting-multidimension
clc; clear; 
para.dt=1e-1;
opt.isfixN=0; % ==0: stop until converge
              % ==1: stop after a fixed number of steps
% parameter setup (default)
if ~isfield(para,'NN'), para.NN=0; end; % # of grid points per dimension
if ~isfield(para,'R'), para.R=16; end;
if ~isfield(para,'beta'), para.beta=10; end;
if ~isfield(para,'delta'), para.delta=10; end;

opt.NNmax=4;
NN=para.NN;

while(para.NN<=opt.NNmax+1)
    
    para.dh=1/2^para.NN; para.M=2*para.R/para.dh; 
    para.x=linspace(-para.R,para.R,para.M+1);
    para.y=para.x; 
    para.X=kron(para.x(2:para.M),ones(para.M-1,1));
    para.Y=kron(ones(para.M-1,1),para.y(2:para.M))';
    para.V=(para.X.^2+para.Y.^2)/2;
    para.v=mat2vec_2D(para.V,para.M-1);
    
    I = speye(para.M-1);
    e = ones(para.M-1,1);
    D2 = spdiags([e -2*e e]/para.dh^2, [-1 0 1], para.M-1, para.M-1);
    para.Delta = kron(D2,I)+kron(I,D2);

    % local v.r.
    x=para.x; y=para.y; 
    v=para.v; M=para.M; dh=para.dh;  
    beta=para.beta; delta=para.delta;
    % initial value setup
    if (para.NN==NN)
        gm=2;
        PHI=sqrt(gm)/pi^(1/2)*exp(-gm*para.X.*para.X/2-gm*para.Y.*para.Y/2);
    else
        PHI0=finer_2D(PHI0,para);
        PHI=PHI0(2:end-1,2:end-1);
    end
    phi=mat2vec_2D(PHI,para.M-1);
    % get ground state
    [phi,~]=ground_solver_2D(phi,para,opt);
    PHI=reshape(phi,M-1,M-1)'; % a=mat2vec_2D(A,M); A=reshape(a,M,M)';
    PHI0=zeros(M+1,M+1);
    PHI0(2:M,2:M)=PHI;
    dPHIx=(PHI0(2:M+1,:)-PHI0(1:M,:))/dh;
    dPHIy=(PHI0(:,2:M+1)-PHI0(:,1:M))/dh;
    E=dh^2*(0.5*sum(sum(dPHIx.^2))+0.5*sum(sum(dPHIy.^2))+sum(v.*phi.^2)+0.5*beta*sum(phi.^4)+...
        2*delta*sum(sum(PHI0(1:M,:).^2.*dPHIx.^2))+2*delta*sum(sum(PHI0(:,1:M).^2.*dPHIy.^2)));
    disp(E)
  
    para.NN=para.NN+1;
end
para.NN=para.NN-1;

% save .mat
fname = strcat('GPE-FD2d-Bet-',num2str(beta),'-Del-',num2str(delta),'-NN-',num2str(para.NN));
save(strcat(fname,'.mat'),'E','PHI','x','y','para');



figure
image(x,y,PHI,'CDataMapping','scaled')
colormap(jet)
axis off;
view(2)

