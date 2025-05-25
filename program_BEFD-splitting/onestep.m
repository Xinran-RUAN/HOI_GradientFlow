function[phi,M_err,E_list]=onestep(phi,para,E_list)
    
    beta=para.beta; delta=para.delta;
    M=para.M; dh=para.dh; dt=para.dt;
    V=para.V;

    temp1=phi;
    % construct dphi
    dphi=(temp1(3:M+1)-temp1(1:M-1))/(2*dh);

    % construct the vector a,b,c,d
    const_s=(0.5+2*delta*temp1(2:M).^2)/((dh^2))*dt;
    a=-const_s; c=-const_s;
    % explicit 
    if beta>0
        b=1+const_s*2+dt*V(2:M)+dt*beta*temp1(2:M).^2; 
        d=temp1(2:M)+dt*2*delta*dphi.^2.*temp1(2:M);
    else
        b=1+const_s*2+dt*V(2:M); 
        d=temp1(2:M)+dt*2*delta*dphi.^2.*temp1(2:M)-dt*beta*temp1(2:M).^2.*temp1(2:M);
    end

%%
    %call TDMA_solver to get updated solution
    temp2=TDMA(a,b,c,d,M-1);
    % Remark: use iterative methods for general (multi-D) problems

%%  
    %normalize temp2
    l2_norm=sum(temp2.^2*dh);
    
    temp2=temp2/sqrt(l2_norm);
    %update temp1 and compute the error
    M_err=sqrt(sum((temp1(2:M)-temp2).^2)*dh);
    % test energy
    phi=[0;temp2;0];
    dphi=(phi(3:M+1)-phi(1:M-1))/(2*dh);
    
    E=dh*(sum(0.5*dphi.^2)+sum(V.*phi.^2)+0.5*beta*sum(phi.^4)+sum(2*delta*phi(2:M).^2.*dphi.^2));
    disp(M_err)
    
    E_list=[E_list,E];