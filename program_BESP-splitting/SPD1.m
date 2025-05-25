% compute the first derivative
% INPUT: M=length(f), R=length(interval)
function [Df]=SPD1(f,M,R)
xi=(-M/2:M/2-1)'*pi/R;
Ff=myfft(f);
FDf=1i*Ff.*xi;
Df=myifft(FDf);
