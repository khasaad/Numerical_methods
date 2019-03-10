%%% Traitement explicite du terme de r?action %%%
clear all; clc
%% D?finir les variables
xmin=0;              %%minimum valeur de x
xmax=4;              %% L = 4, maximum valeur de x. Longueur du domaine
M=200;               %%Numero de pas d'espace 
dt=0.01;             %%pas de temps
t=0;                 %%temps
tfinale=10;          %%temps finale
omega=0.9;           %%porosit?
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


%%faire une boucle en temps
while(t<=tfinale)
    %%calculer conditions aux limites
    cnp1(1)=t*exp(-t/2);     %%conditions aux limites gauche
    cnp1(M)=0;               %%conditions aux limites droite
    for j = 2 : M -1
        cnp1(j) = c(j) - u*dt/(omega*dx)*(c(j)-c(j-1)) + v*dt/...
        (omega*dx*dx)*(c(j+1) -2*c(j) +c(j-1)) -k*dt/omega*c(j)*s(j);        
    end
    
    for j = 1 : M
        snp1(j) = s(j)-dt*p*k*c(j)*s(j);
    end
    
    %%mettre ? jour t, c et s
    t = t + dt;
    fprintf('temps de simulation       = %f \n ',t);
    c = cnp1;    
    s = snp1;    
end

figure(2)
colstr = {'r--', '-g', 'm', 'k--'};
plot(x,c,colstr{1}, 'linewidth', 2)
hold on;
plot(x,s,colstr{2}, 'linewidth', 2)
title('explicite', 'fontsize', 18)
legend('c','s')
xlabel('x', 'fontsize', 18)
ylabel('c', 'fontsize', 18)
set(gca, 'fontsize', 18)
    
