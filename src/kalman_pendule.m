%% kalman_pendule.m — LQR avec mesures bruitées : dérivation vs filtre de Kalman (LQG)
% Seules x et theta sont mesurées (codeurs), avec du bruit. Les vitesses doivent
% être estimées. Trois versions de la même commande LQR sont comparées :
%   1) état réel complet (cas idéal, référence)
%   2) vitesses obtenues en dérivant les mesures bruitées
%   3) état estimé par un filtre de Kalman (LQG)
% Commande numérique échantillonnée à 100 Hz, moteur saturé à 10 N,
% et une poussée de 5 N pendant 0,1 s à t = 2 s.
%
% Prérequis : derive_eom.m (cartpole_dyn.m, cartpole_lin.m).
% Requiert  : Control System Toolbox.

clear; clc; close all;

%% 1) Modèle et LQR (identiques à lqr_pendule.m)
M = 0.5;  m = 0.2;  l = 0.3;  J = m*(2*l)^2/12;  b = 0.1;  g = 9.81;
p = [M m l J b g];
[A, B] = cartpole_lin(p);
K = lqr(A, B, diag([10 100 1 1]), 1);

%% 2) Échantillonnage et capteurs
Ts   = 0.01;                    % période de la commande : 100 Hz
dt   = 1e-3;                    % pas de simulation du système réel
nsub = round(Ts/dt);
C    = [1 0 0 0;                % on mesure x...
        0 1 0 0];               % ...et theta, pas les vitesses

sig_x  = 0.5e-3;                % bruit du codeur du chariot : 0,5 mm
sig_th = deg2rad(0.2);          % bruit du codeur du pendule : 0,2°

%% 3) Filtre de Kalman (régime permanent)
sysd = c2d(ss(A, B, eye(4), 0), Ts);           % modèle discret
Ad = sysd.A;  Bd = sysd.B;

Rv    = diag([sig_x^2, sig_th^2]);              % confiance dans les mesures
sig_w = 2;                                    % perturbations de force inconnues [N]
Qw    = Bd*Bd.'*sig_w^2 + diag([1e-9 1e-9 1e-7 1e-7]);   % confiance dans le modèle

try
    P = idare(Ad.', C.', Qw, Rv);               % covariance d'erreur en régime permanent
catch
    P = dare(Ad.', C.', Qw, Rv);                % versions plus anciennes de MATLAB
end
Lk = P*C.' / (C*P*C.' + Rv);                    % gain de Kalman

pe = log(eig(Ad*(eye(4) - Lk*C))) / Ts;         % pôles de l'estimateur (équivalent continu)
disp('Gain de Kalman Lk :'); disp(Lk)
disp('Pôles de l''estimateur :'); disp(pe)
disp('Pôles de la commande  :'); disp(eig(A - B*K))

%% 4) Simulation des trois versions
Fmax = 10;  T = 6;  N = round(T/Ts);
X0   = [0; deg2rad(10); 0; 0];
pousse = @(t) 5*(t >= 2 & t < 2.1);            % coup de 5 N pendant 0,1 s

noms = {'état réel (idéal)', 'dérivation des mesures', 'filtre de Kalman'};
t    = (0:N-1).' * Ts;
Xv   = zeros(N, 4, 3);                          % état réel
Xe   = zeros(N, 4, 3);                          % état utilisé par la commande
Fv   = zeros(N, 3);

for mode = 1:3
    rng(1);                                     % même bruit pour les trois versions
    X  = X0;
    xh = [C*X0; 0; 0];                          % Kalman : positions mesurées, vitesses nulles
    yprec = C*X0;
    u  = 0;

    for k = 1:N
        y = C*X + [sig_x; sig_th].*randn(2,1);  % mesure bruitée

        switch mode
            case 1                              % état réel
                xe = X;
            case 2                              % dérivation numérique
                xe = [y; (y - yprec)/Ts];
            case 3                              % Kalman : prédiction puis correction
                xm = Ad*xh + Bd*u;
                xh = xm + Lk*(y - C*xm);
                xe = xh;
        end
        yprec = y;

        u = max(min(-K*xe, Fmax), -Fmax);       % commande, maintenue pendant Ts
        Xv(k,:,mode) = X.';  Xe(k,:,mode) = xe.';  Fv(k,mode) = u;

        for s = 1:nsub                          % système réel entre deux échantillons
            ts = t(k) + (s-1)*dt;
            F  = u + pousse(ts);
            k1 = cartpole_dyn(X,           F, p);
            k2 = cartpole_dyn(X + dt/2*k1, F, p);
            k3 = cartpole_dyn(X + dt/2*k2, F, p);
            k4 = cartpole_dyn(X + dt*k3,   F, p);
            X  = X + dt/6*(k1 + 2*k2 + 2*k3 + k4);
        end
    end
end

%% 5) Indicateurs
calme  = t >= 4;                                % après le retour au calme
apres  = t >= 2 & t < 3.5;                      % réaction à la poussée
rmsv   = @(v) sqrt(mean(v.^2));
fprintf('\n%-24s | theta max après poussée | bruit sur theta | bruit sur F | erreur sur dtheta\n', '');
for mode = 1:3
    th  = rad2deg(Xv(:,2,mode));
    edt = Xe(:,4,mode) - Xv(:,4,mode);
    fprintf('%-24s |        %5.2f°          |     %5.3f°     |   %5.2f N   |   %5.3f rad/s\n', ...
        noms{mode}, max(abs(th(apres))), std(th(calme)), std(Fv(calme,mode)), rmsv(edt(t >= 1)));
end

%% 6) Tracés
coul = [0 0 0; 0.85 0.33 0.1; 0 0.45 0.74];

figure('Color', 'w', 'Name', 'Comparaison', 'Position', [100 100 820 650]);
subplot(2,1,1); hold on; grid on;
for mode = [2 3 1], plot(t, rad2deg(Xv(:,2,mode)), 'Color', coul(mode,:), 'LineWidth', 1.1); end
ylabel('\theta [deg]'); legend(noms([2 3 1]), 'Location', 'northeast');
title('LQR avec mesures bruitées de x et \theta (100 Hz, poussée de 5 N à t = 2 s)');
subplot(2,1,2); hold on; grid on;
for mode = [2 3 1], plot(t, Fv(:,mode), 'Color', coul(mode,:), 'LineWidth', 1.1); end
ylabel('F [N]'); xlabel('t [s]');

figure('Color', 'w', 'Name', 'Vitesse estimée', 'Position', [150 150 820 400]);
hold on; grid on;
plot(t, Xe(:,4,2), 'Color', coul(2,:));
plot(t, Xe(:,4,3), 'Color', coul(3,:), 'LineWidth', 1.4);
plot(t, Xv(:,4,3), 'k--', 'LineWidth', 1.2);
xlim([1.5 3.5]); ylabel('d\theta/dt [rad/s]'); xlabel('t [s]');
legend('dérivation des mesures', 'filtre de Kalman', 'vitesse réelle', 'Location', 'best');
title('Estimation de la vitesse angulaire autour de la poussée');
