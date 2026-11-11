%% lqr_pendule.m — Stabilisation du pendule inversé par LQR
% Étape 2 du projet. Prérequis : avoir lancé derive_eom.m
% (qui crée cartpole_dyn.m et cartpole_lin.m dans le même dossier).
% Requiert : Control System Toolbox (fonctions ctrb et lqr)

clear; clc; close all;

%% 1) Paramètres (les mêmes que dans derive_eom.m)
M = 0.5;  m = 0.2;  l = 0.3;
J = m*(2*l)^2/12;
b = 0.1;  g = 9.81;
p = [M m l J b g];

%% 2) Modèle linéarisé
[A, B] = cartpole_lin(p);

%% 3) Commandabilité : peut-on piloter les 4 états avec la seule force F ?
Co = ctrb(A, B);
fprintf('Rang de commandabilité : %d / 4\n', rank(Co));

%% 4) Réglage LQR : les priorités
%        x    theta   dx   dtheta
Q = diag([10,  100,    1,    1]);   % importance de chaque écart
R = 1;                              % prix de la force
Fmax = 10;   % limite du moteur [N]
K = lqr(A, B, Q, R);
disp('Gain K = [k_x  k_theta  k_dx  k_dtheta] :'); disp(K)
disp('Pôles en boucle fermée :'); disp(eig(A - B*K))

%% 5) Simulation sur le modèle NON linéaire (le "vrai" pendule)
angle0 = 30;                          % angle de départ [deg]
X0 = [0; deg2rad(angle0); 0; 0];        % départ : pendule penché de 10°, tout immobile
loi = @(X) max(min(-K*X, Fmax), -Fmax);                    % la loi de commande F = -K*X

[t, Xs] = ode45(@(t,X) cartpole_dyn(X, loi(X), p), [0 5], X0);
Fs = max(min(-Xs*K.', Fmax), -Fmax);                       % force appliquée à chaque instant

%% 6) Tracés
figure('Name', 'LQR');
subplot(3,1,1); plot(t, Xs(:,1));          grid on; ylabel('x [m]');
title(sprintf('Stabilisation LQR (départ à %d°, Fmax = %d N)', angle0, Fmax));subplot(3,1,2); plot(t, rad2deg(Xs(:,2))); grid on; ylabel('\theta [deg]');
subplot(3,1,3); plot(t, Fs);               grid on; ylabel('F [N]'); xlabel('t [s]');

fprintf('Force maximale : %.1f N\n', max(abs(Fs)));
