% only for n to be even
function [Ff]=myfft(f)
[m,n]=size(f);
% Fourier transform
if (m-1)*(n-1)~=0
    disp('Error')
end
if(m==1)
    g1=fft(f)/n;
    Ff=[g1(n/2+1:n),g1(1:n/2)];
else
    g1=fft(f)/m;
    Ff=[g1(m/2+1:m);g1(1:m/2)];
end
end