%% lqr_double.m — Stabilisation du double pendule inversé par LQR
% Étape 3 du projet. Prérequis : avoir lancé derive_eom_nlink.m avec N = 2
% (qui crée nlink_M.m, nlink_h.m, nlink_linM.m et nlink_model.mat).
% Requiert : Control System Toolbox (fonctions ctrb et lqr)

clear; clc; close all;

%% 1) Récupérer le modèle généré par derive_eom_nlink.m
load('nlink_model.mat', 'A', 'B', 'p', 'N');
nq = N + 1;                                  % nombre de degrés de liberté (chariot + N barres)
fprintf('Modèle chargé : N = %d barres, %d états\n', N, 2*nq);

%% 2) Commandabilité
fprintf('Rang de commandabilité : %d / %d\n', rank(ctrb(A, B)), 2*nq);

%% 3) Réglage LQR
%        x     theta1..thetaN       vitesses
Q = diag([10,  100*ones(1, N),      ones(1, nq)]);
R = 1;
Fmax = 10;                                  % pas de limite moteur pour commencer

K = lqr(A, B, Q, R);
disp('Gain K = [k_x  k_th1  k_th2  k_dx  k_dth1  k_dth2] :'); disp(K)
disp('Pôles en boucle fermée :'); disp(eig(A - B*K))

%% 4) Modèle non linéaire complet :  M(q)*qdd = Bu*F - D*qd - h(q,qd)
bc  = p(2);                                  % frottement du chariot
dyn = @(X, F) [ X(nq+1:end);
                nlink_M(X(1:nq), p) \ ( [F; zeros(N,1)] ...
                                      - [bc*X(nq+1); zeros(N,1)] ...
                                      - nlink_h(X(1:nq), X(nq+1:end), p) ) ];

loi = @(X) max(min(-K*X, Fmax), -Fmax);      % F = -K*X, saturée à +/- Fmax

%% 5) Simulation
angle0 = 20;                                  % inclinaison de départ des deux barres [deg]
X0 = zeros(2*nq, 1);
X0(2:nq) = deg2rad(angle0);                  % tout le double pendule penché d'un bloc

[t, Xs] = ode45(@(t,X) dyn(X, loi(X)), [0 5], X0);
Fs = max(min(-Xs*K.', Fmax), -Fmax);

%% 6) Tracés
figure('Name', 'LQR double pendule');
subplot(3,1,1); plot(t, Xs(:,1)); grid on; ylabel('x [m]');
title(sprintf('Double pendule : LQR (départ à %d°, Fmax = %g N)', angle0, Fmax));
subplot(3,1,2); plot(t, rad2deg(Xs(:,2:nq))); grid on; ylabel('\theta [deg]');
legend(compose('\\theta_%d', 1:N));
subplot(3,1,3); plot(t, Fs); grid on; ylabel('F [N]'); xlabel('t [s]');

fprintf('Force maximale : %.1f N\n', max(abs(Fs)));
fprintf('Déplacement max du chariot : %.2f m\n', max(abs(Xs(:,1))));
