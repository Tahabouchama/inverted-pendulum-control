%% traj_opt_double.m — Swing-up du double pendule : optimisation de trajectoire
% Cherche la force F(t) qui amène les deux barres de la position basse à la
% position haute en T secondes, avec le moins d'effort possible.
% Méthode : collocation directe (règle des trapèzes) résolue par fmincon.
%
% Prérequis : avoir lancé derive_eom_nlink.m avec N = 2
%             (crée nlink_M.m, nlink_h.m et nlink_model.mat).
% Requiert  : Optimization Toolbox (fmincon).
% Durée     : quelques minutes.

clear; clc; close all;

if ~license('test', 'Optimization_Toolbox')
    error('Optimization Toolbox non disponible : fmincon est nécessaire.');
end

%% 1) Modèle
S  = load('nlink_model.mat', 'p', 'N');
p  = S.p;  N = S.N;
nq = N + 1;  ns = 2*nq;                       % 3 degrés de liberté, 6 états

%% 2) Réglages du problème
T    = 3.0;      % durée de la manœuvre [s]
Nn   = 61;       % nombre d'instants (noeuds)
Fmax = 20;       % force maximale [N]
xmax = 1.0;      % demi-longueur des rails [m]
h    = T/(Nn-1); % pas de temps entre deux noeuds
t    = linspace(0, T, Nn);

X0 = [0; pi; pi; 0; 0; 0];                    % départ : les deux barres pendues
Xf = [0;  0;  0; 0; 0; 0];                    % arrivée : les deux barres debout

%% 3) Bornes : limites physiques + départ et arrivée imposés
LBX = -inf(ns, Nn);   UBX = inf(ns, Nn);
LBX(1,:) = -xmax;     UBX(1,:) = xmax;        % rails
LBX(:,1)   = X0;      UBX(:,1)   = X0;        % état initial fixé
LBX(:,end) = Xf;      UBX(:,end) = Xf;        % état final fixé
lb = [LBX(:); -Fmax*ones(Nn,1)];
ub = [UBX(:);  Fmax*ones(Nn,1)];

%% 4) Point de départ de l'optimiseur : interpolation linéaire bas -> haut
s   = linspace(0, 1, Nn);
Xg  = X0 + (Xf - X0)*s;
Xg(nq+2:end, :) = -pi/T;                      % vitesses angulaires cohérentes
z0  = [Xg(:); zeros(Nn,1)];

%% 5) Résolution
opts = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'SpecifyObjectiveGradient', true, ...
    'MaxIterations', 500, ...
    'MaxFunctionEvaluations', 5e5, ...
    'ConstraintTolerance', 1e-6, ...
    'Display', 'iter');

fprintf('Optimisation en cours (quelques minutes)...\n');
tic;
[z, Jopt, exitflag] = fmincon(@(z) cout(z, ns, Nn, h), z0, [], [], [], [], lb, ub, ...
                              @(z) physique(z, ns, Nn, h, p, N), opts);
fprintf('\nTerminé en %.0f s  |  exitflag = %d  |  effort = %.1f N²·s\n', toc, exitflag, Jopt);
if exitflag <= 0
    fprintf('Aucune trajectoire valide trouvée : essayer T = 3.5 ou Fmax = 30.\n');
end

Xopt = reshape(z(1:ns*Nn), ns, Nn);
Fopt = z(ns*Nn+1:end).';
fprintf('Force max : %.1f N  |  Déplacement max du chariot : %.2f m\n', ...
        max(abs(Fopt)), max(abs(Xopt(1,:))));

save('traj_double.mat', 't', 'Xopt', 'Fopt', 'p', 'N', 'T');

%% 6) Rejeu en boucle OUVERTE : on applique F(t) tel quel au modèle complet
dt = 1e-3;  tr = 0:dt:T;  nr = numel(tr);
Xr = zeros(ns, nr);  Xr(:,1) = X0;
for k = 1:nr-1
    Fk = interp1(t, Fopt, tr(k));
    Fm = interp1(t, Fopt, tr(k) + dt/2);
    Fn = interp1(t, Fopt, tr(k+1));
    X  = Xr(:,k);
    k1 = dyn(X,           Fk, p, N);
    k2 = dyn(X + dt/2*k1, Fm, p, N);
    k3 = dyn(X + dt/2*k2, Fm, p, N);
    k4 = dyn(X + dt*k3,   Fn, p, N);
    Xr(:,k+1) = X + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end
fprintf('Rejeu boucle ouverte, écart final : θ1 = %+.1f°, θ2 = %+.1f°\n', ...
        rad2deg(Xr(2,end)), rad2deg(Xr(3,end)));

%% 7) Tracés
figure('Color', 'w', 'Name', 'Trajectoire optimale', 'Position', [100 100 820 650]);
subplot(3,1,1);
plot(t, rad2deg(Xopt(2,:)), 'o-', t, rad2deg(Xopt(3,:)), 'o-', 'MarkerSize', 3); hold on;
plot(tr, rad2deg(Xr(2,:)), '--', tr, rad2deg(Xr(3,:)), '--');
grid on; ylabel('\theta [deg]');
legend('\theta_1 planifié', '\theta_2 planifié', '\theta_1 rejoué', '\theta_2 rejoué', ...
       'Location', 'best');
title('Swing-up du double pendule : trajectoire optimale');
subplot(3,1,2); plot(t, Xopt(1,:), 'o-', 'MarkerSize', 3); grid on; ylabel('x [m]');
subplot(3,1,3); plot(t, Fopt, 'o-', 'MarkerSize', 3); grid on; ylabel('F [N]'); xlabel('t [s]');

%% ---------- Fonctions locales (tout passe par les arguments) ----------
function [J, g] = cout(z, ns, Nn, h)
% Effort : somme des F² (règle des trapèzes), et son gradient exact.
    F = z(ns*Nn+1:end);
    w = h*ones(Nn,1);  w([1 end]) = h/2;
    J = sum(w .* F.^2);
    g = [zeros(ns*Nn,1); 2*w.*F];
end

function [c, ceq] = physique(z, ns, Nn, h, p, N)
% Contraintes de dynamique : entre deux noeuds, l'état doit suivre le modèle
% (règle des trapèzes). ceq = 0 quand la physique est respectée.
    X = reshape(z(1:ns*Nn), ns, Nn);
    F = z(ns*Nn+1:end);
    fx = zeros(ns, Nn);
    for k = 1:Nn
        fx(:,k) = dyn(X(:,k), F(k), p, N);
    end
    def = X(:,2:end) - X(:,1:end-1) - h/2*(fx(:,2:end) + fx(:,1:end-1));
    ceq = def(:);
    c   = [];
end

function dX = dyn(X, F, p, N)
% Modèle non linéaire complet :  M(q)*qdd = [F - b*dx; 0...] - h(q,qd)
    nq = N + 1;
    q  = X(1:nq);  dq = X(nq+1:end);
    Q  = [F - p(2)*dq(1); zeros(N,1)];
    dX = [dq; nlink_M(q, p) \ (Q - nlink_h(q, dq, p))];
end
