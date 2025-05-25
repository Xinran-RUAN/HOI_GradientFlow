% solve for one step: \phi^n->\phi^{n+1}
% via BiCGSTAB method 
% INPUT: phi_o=\phi^n, dphi=dphi^n
% OUTPUT: phi_n=\phi^{n+1}
function[phi_n,M_err,E_list]=onestep(phi_o,dphi,para,E_list)
% get parameters from para
beta=para.beta; delta=para.delta;
M=para.M;  V=para.V; 
dt=para.dt; dh=para.dh; 
R=para.R;  xi=(-M/2:M/2-1)'*pi/R;

% set initial value and parameters
x0=phi_o; 
b=phi_o/dt+2*delta*abs(dphi).* abs(dphi).*phi_o;
absphi2=abs(phi_o).*abs(phi_o);
Ax0=x0/dt+V.*x0+beta*absphi2.*x0+(0.5+2*delta*absphi2).*myifft(myfft(x0).*xi.^2); % temp=A*x0
r=b-Ax0; sr0=r; % sr0: star_r0
p0=r;

err=1;

while (err>1e-8)
    % update from \phi^{+,m} to \phi^{+,m+1}  
    % update x0, r, p0
    %=========================================================
    %   BiCGSTAB.  
    %=========================================================
    Ap0=p0/dt+V.*p0+beta*absphi2.*p0+(0.5+2*delta*absphi2).*myifft(myfft(p0).*xi.^2); % temp=A*p0
    temp=sum(r.*sr0);
    alpha=temp/sum(Ap0.*sr0);
    s=r-alpha*Ap0;
    As=s/dt+V.*s+beta*absphi2.*s+(0.5+2*delta*absphi2).*myifft(myfft(s).*xi.^2);
    omega=sum(As.*s)/sum(As.*As);
    xn=x0+alpha*p0+omega*s;
    rn=s-omega*As;
    bt=sum(rn.*sr0)/temp*alpha/omega;
    pn=rn+bt*(p0-omega*Ap0);
    
    err=max(abs(rn));
    p0=pn; x0=xn; r=rn;
end





% normalize step
l2_norm=sum(x0.^2)*dh;
x0=x0/sqrt(l2_norm);
%update temp1 and compute the error
M_err=sqrt(sum((abs(x0-phi_o)).^2)*dh);
phi_n=x0;



% ! E: spectral
phi=phi_n; rho=phi.^2;
fphi=myfft(phi);
frho=myfft(rho);
E=R*sum(abs(fphi).^2.*xi.^2)+delta*R*sum(abs(frho).^2.*xi.^2)+dh*(sum(V.*phi.^2)+0.5*beta*sum(phi.^4));
% Note: dh*sum(0.5*dphi.^2)=R*sum(fphi.^2.*xi.^2); 
E_list=[E_list,E];
