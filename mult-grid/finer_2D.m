% Linear interpolation
function[PHI]=finer_2D(PHI,para)
M=2*para.R*2^para.NN; 
temp=zeros(M+1,M+1);
temp(1:2:end,1:2:end)=PHI;
temp(2:2:end-1,1:2:end)=0.5*(PHI(1:end-1,:)+PHI(2:end,:));
temp(1:2:end,2:2:end-1)=0.5*(PHI(:,1:end-1)+PHI(:,2:end));
temp(2:2:end-1,2:2:end-1)=0.25*(PHI(1:end-1,1:end-1)+PHI(2:end,1:end-1)+...
    PHI(1:end-1,2:end)+PHI(2:end,2:end));
PHI=temp;