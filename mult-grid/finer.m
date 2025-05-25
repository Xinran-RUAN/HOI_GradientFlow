% Linear interpolation
function[temp]=finer(phi,para)
M=2*para.R*2^para.NN; 
temp=zeros(M+1,1);
temp(1:2:end)=phi;
temp(2:2:end-1)=0.5*(phi(1:end-1)+phi(2:end));
