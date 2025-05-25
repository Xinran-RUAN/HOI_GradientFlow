% This file is used for getting the ground state
% for beta>=0
% Input: M+1=length(phi),
function[phi,E_list]=ground_solver_2D(phi,para,opt)
E_list=[]; 
if opt.isfixN==1
    for jj=1:opt.Nmax
        [phi,~]=onestep_2D(phi,para);
    end
else
    M_err=1;  vep=1e-6;
    while(M_err/para.dt>vep)
        [phi,M_err]=onestep_2D(phi,para);
    end
end

