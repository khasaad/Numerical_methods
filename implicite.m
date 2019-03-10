%%d?finir les variables
xmin=0;              %%minimum valeur de x
xmax=4;              %%maximum valeur de x. Longueur du domaine
M=200;               %%Numero de pas d'espace 
dt=0.01;             %%pas de temps
t=0;                 %%temps
tfinale=10;          %%temps finale
omega=0.9;           %%porosit??
u=0.4;               %%vitesse moyenne
v=0.01;              %%coefficient de diffusion
k=5;                 %%constante de reaction
p=10;                %%valence de l'anticorps

%%Discretisation
dx = (xmax - xmin)/(M-1);  %%pas de space
x=xmin : dx : xmax;        %%pr?cisant les points de la grille (maillage)

%%conditions initiales
c(1:M)=0;
s(1:M)=10;

while(t<=tfinale)
    
    sj = s(2:M-1);
    cj = c(2:M-1);
    cjm1=c(1:M-2);
    cjp1=c(3:M);
    
    fhandle = @fonction_anticorps;
    options = optimset('display','off');
    cons = fsolve(fhandle,cj,options,omega,dt,cj,u,dx,cjm1,v,cjp1,k,sj,p);
    cnp1 = [t*exp(-t/2),cons,0];
    snp1 = s./(1+dt*p*k*cnp1);
    %% mettre ? jour t, c et s
    t = t + dt;
    fprintf('temps de simulation       = %f  \n',t);
    c = cnp1;    
    s = snp1;
end

figure(1)
colstr = {'r--', '-g', 'm', 'k--'};
plot(x,c,colstr{1}, 'linewidth', 2)
hold on;
plot(x,s,colstr{2}, 'linewidth', 2)
title('implicite', 'fontsize', 18)
legend('c','s')
xlabel('x', 'fontsize', 18)
ylabel('c', 'fontsize', 18)
set(gca, 'fontsize', 18)
    
