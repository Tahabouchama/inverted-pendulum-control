%% derive_eom.m — Pendule inversé sur chariot : modèle de Lagrange
% Étape 1 du projet portfolio.
% Dérive les équations du mouvement, génère la dynamique non linéaire
% et le modèle linéarisé autour de l'équilibre haut (theta = 0).
%
% Conventions :
%   x     : position du chariot [m]
%   theta : angle du pendule depuis la verticale HAUTE [rad] (theta=0 : debout)
%   F     : force appliquée au chariot [N]
%   Etat  : X = [x; theta; dx; dtheta]
%
% Requiert : Symbolic Math Toolbox

clear; clc;

%% 1) Variables symboliques
syms x th dx dth ddx ddth F real
syms M m l J b g real    % M chariot, m pendule, l dist. pivot->CdG, J inertie au CdG, b frottement chariot

%% 2) Énergies
% Position du centre de gravité du pendule
xp  = x + l*sin(th);
yp  = l*cos(th);
dxp = dx + l*cos(th)*dth;
dyp = -l*sin(th)*dth;

T = 1/2*M*dx^2 + 1/2*m*(dxp^2 + dyp^2) + 1/2*J*dth^2;   % cinétique
V = m*g*yp;                                             % potentielle
L = T - V;                                              % lagrangien

%% 3) Équations de Lagrange : d/dt(dL/dq') - dL/dq = Q
q   = [x; th];
dq  = [dx; dth];
ddq = [ddx; ddth];

dLdq      = jacobian(L, q).';
dLddq     = jacobian(L, dq).';
ddt_dLddq = jacobian(dLddq, [q; dq]) * [dq; ddq];   % dérivée totale par rapport au temps

Q   = [F - b*dx; 0];                 % forces généralisées (commande + frottement visqueux)
eqs = simplify(ddt_dLddq - dLdq - Q);

disp('Équations du mouvement (= 0) :'); pretty(eqs)

%% 4) Forme d'état non linéaire  dX/dt = f(X, F)
sol = solve(eqs == 0, ddq);
X   = [x; th; dx; dth];
f   = simplify([dx; dth; sol.ddx; sol.ddth]);

%% 5) Linéarisation autour de l'équilibre haut (X = 0, F = 0)
A = simplify(subs(jacobian(f, X), [X; F], zeros(5,1)));
B = simplify(subs(jacobian(f, F), [X; F], zeros(5,1)));
disp('A ='); pretty(A)
disp('B ='); pretty(B)

%% 6) Génération des fonctions MATLAB réutilisables
p_sym = [M m l J b g];
matlabFunction(f,    'File', 'cartpole_dyn', 'Vars', {X, F, p_sym});
matlabFunction(A, B, 'File', 'cartpole_lin', 'Vars', {p_sym});
fprintf('\nFichiers générés : cartpole_dyn.m, cartpole_lin.m\n');

%% 7) Vérifications numériques
% Paramètres (ordre de grandeur d'un banc type Quanser)
Mv = 0.5;  mv = 0.2;  lv = 0.3;      % tige de longueur 2l, CdG au milieu
Jv = mv*(2*lv)^2/12;                 % inertie d'une tige au CdG
bv = 0.1;  gv = 9.81;
p  = [Mv mv lv Jv bv gv];

[An, Bn] = cartpole_lin(p);
disp('Pôles en boucle ouverte :'); disp(eig(An))
% Attendu : un pôle réel positif (instabilité), un pôle nul (position libre)

% Simulation non linéaire en boucle ouverte, petit écart initial
X0 = [0; 0.05; 0; 0];
[t, Xs] = ode45(@(t,X) cartpole_dyn(X, 0, p), [0 3], X0);

figure('Name','Boucle ouverte');
subplot(2,1,1); plot(t, Xs(:,1)); grid on; ylabel('x [m]');
title('Boucle ouverte : le pendule tombe (attendu)');
subplot(2,1,2); plot(t, rad2deg(Xs(:,2))); grid on; ylabel('\theta [deg]'); xlabel('t [s]');

% Contrôle de cohérence : sans frottement, l'énergie totale doit se conserver
p0 = p; p0(5) = 0;
[t0, X0s] = ode45(@(t,X) cartpole_dyn(X, 0, p0), [0 3], X0, odeset('RelTol',1e-9,'AbsTol',1e-9));
E = 1/2*(Mv+mv)*X0s(:,3).^2 + mv*lv*X0s(:,3).*X0s(:,4).*cos(X0s(:,2)) ...
  + 1/2*(Jv+mv*lv^2)*X0s(:,4).^2 + mv*gv*lv*cos(X0s(:,2));
fprintf('Dérive d''énergie max (b=0) : %.2e J\n', max(abs(E - E(1))));
