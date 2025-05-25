% TDMA: tridiagonal matrix algorithm
% a_i*x_{i-1}+b_i*x_i+c_i*x_{i+1}=d_i, i=1,2,...N
% a(1)=0, c(N)=0
function[x]=TDMA(a,b,c,d,N)
x=zeros(size(d));
if N==1
    x=d/b;
else
    % forward sweeping
    a(1)=0; c(N)=0;
    c(1)=c(1)/b(1);
    d(1)=d(1)/b(1);
    for ii=2:N-1
        temp=b(ii)-c(ii-1)*a(ii);
        c(ii)=c(ii)/temp;
        d(ii)=(d(ii)-d(ii-1)*a(ii))/temp;
    end
    d(N)=(d(N)-d(N-1)*a(N))/(b(N)-c(N-1)*a(N));
    % backward sweeping
    x(N)=d(N);
    for ii=N-1:-1:1
        x(ii)=d(ii)-c(ii)*x(ii+1);
    end
end