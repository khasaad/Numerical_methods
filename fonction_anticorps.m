%% fonction non linear pour calculer la concentration d'anticorps
function [g] = fonction_anticorps(cnp1,omega,dt,cj,u,dx,cjm1,v,cjp1,k,sj,p)
    g = omega/dt*(cnp1-cj) + u/dx*(cj-cjm1) - v/(dx*dx)*(cjp1 -2*cj +cjm1) + k*cnp1.*sj./(1+dt*p*k*cnp1);
end

