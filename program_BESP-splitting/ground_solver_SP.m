% This file is used for getting the ground state
% for beta>=0
function[phi_n,E_list]=ground_solver_SP(phi,para,opt)
R=para.R; M=para.M; dt=para.dt;
E_list=[];
%=============iPiano acceleration===============
bt_iPiano=0; % no acceleration when bt=0
phi_o=phi;
%===============================================

if opt.isfixN==1
    for jj=1:opt.Nmax
        % construct dphi
        dphi=SPD1(phi,M,R);
        [phi_n,M_err,E_list]=onestep(phi,dphi,para,E_list);
        phi=phi_n;
        disp(M_err)
    end
else
    M_err=1;  vep=1e-6;
    while(M_err/dt>vep)
        %==================iPiano acceleration===============
        tphi=phi+bt_iPiano*(phi-phi_o);
        %====================================================
        % construct dphi
        dphi=SPD1(tphi,M,R);
        [phi_n,M_err,E_list]=onestep(tphi,dphi,para,E_list);
        phi_o=phi; phi=phi_n; 
        disp(E_list(end))         
    end
end