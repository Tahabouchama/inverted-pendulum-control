%% tvlqr_double.m — Swing-up du double pendule : suivi de trajectoire par TVLQR
% Étape 2 du swing-up. Suit la trajectoire trouvée par traj_opt_double.m
% avec un retour d'état variant dans le temps, puis bascule sur le LQR fixe
% de la position haute.
%
%   F(t) = F_plan(t) - K(t) * (X(t) - X_plan(t))     pendant la manoeuvre
%   F(t) = -K_haut * X(t)                            ensuite
%
% Prérequis : traj_double.mat (traj_opt_double.m) et nlink_model.mat,
%             nlink_M.m, nlink_h.m (derive_eom_nlink.m avec N = 2).
% Requiert  : Control System Toolbox (lqr).

clear; clc; close all;

%% 1) Chargement de la trajectoire et du modèle
load('traj_double.mat', 't', 'Xopt', 'Fopt', 'p', 'N', 'T');
S  = load('nlink_model.mat', 'A', 'B');
nq = N + 1;  ns = 2*nq;

Xplan = @(tt) interp1(t, Xopt.', tt, 'pchip').';     % état planifié
Fplan = @(tt) interp1(t, Fopt,   tt, 'linear');      % force planifiée

%% 2) Réglages
Q    = diag([10, 100, 100, 1, 1, 1]);   % mêmes priorités que lqr_double.m
R    = 1;
Fmax = 10;                              % moteur d'origine
dth0 = 0;                               % perturbation initiale de theta1 [deg] (essayer 3)
Tfin = T + 2;                           % on observe 2 s après la manoeuvre

%% 3) LQR de la position haute : gain final et condition finale de Riccati
[Khaut, Shaut] = lqr(S.A, S.B, Q, R);

%% 4) Linéarisation le long de la trajectoire
tg  = linspace(0, T, 301);              % grille de calcul des gains
dtg = tg(2) - tg(1);
Ag  = zeros(ns, ns, numel(tg));
Bg  = zeros(ns, numel(tg));
for k = 1:numel(tg)
    [Ag(:,:,k), Bg(:,k)] = linearise(Xplan(tg(k)), Fplan(tg(k)), p, N);
end

%% 5) Riccati à rebours, de t = T jusqu'à t = 0, en partant de Shaut
fprintf('Calcul des gains K(t)...\n');
ric = @(tt, s) riccati(tt, s, tg, dtg, Ag, Bg, Q, R, ns);
[~, Sv] = ode45(ric, fliplr(tg), Shaut(:), odeset('RelTol', 1e-6, 'AbsTol', 1e-8));
Sv = flipud(Sv);                         % remis dans l'ordre chronologique

Kg = zeros(numel(tg), ns);
for k = 1:numel(tg)
    Sk = reshape(Sv(k,:), ns, ns);
    Kg(k,:) = (Bg(:,k).' * Sk) / R;
end
fprintf('Écart entre K(T) et le LQR haut : %.2e (doit être ~0)\n', norm(Kg(end,:) - Khaut));

%% 6) Simulation en boucle fermée (RK4, pas fixe)
dt = 1e-3;  ts = 0:dt:Tfin;  n = numel(ts);
Xs = zeros(ns, n);  Fs = zeros(1, n);  Ms = zeros(1, n);
Xs(:,1) = Xopt(:,1) + [0; deg2rad(dth0); 0; 0; 0; 0];

for k = 1:n
    tk = ts(k);  X = Xs(:,k);
    if tk < T
        Ms(k) = 1;                                       % suivi de trajectoire
        Kk = interp1(tg, Kg, tk);
        F  = Fplan(tk) - Kk*(X - Xplan(tk));
    else
        Ms(k) = 2;                                       % LQR haut
        Xw = X;  Xw(2:nq) = mod(X(2:nq) + pi, 2*pi) - pi;
        F  = -Khaut*Xw;
    end
    Fs(k) = max(min(F, Fmax), -Fmax);

    if k < n
        k1 = dyn(X,           Fs(k), p, N);
        k2 = dyn(X + dt/2*k1, Fs(k), p, N);
        k3 = dyn(X + dt/2*k2, Fs(k), p, N);
        k4 = dyn(X + dt*k3,   Fs(k), p, N);
        Xs(:,k+1) = X + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    end
end

thf = rad2deg(mod(Xs(2:nq,end) + pi, 2*pi) - pi);
fprintf('Angles finaux : theta1 = %+.2f°, theta2 = %+.2f°\n', thf(1), thf(2));
fprintf('Force max : %.1f N  |  Chariot max : %.2f m\n', max(abs(Fs)), max(abs(Xs(1,:))));

%% 7) Tracés
tp = ts(ts <= T);  Xp = Xplan(tp);
figure('Color', 'w', 'Name', 'TVLQR', 'Position', [100 100 820 700]);

subplot(3,1,1); hold on; grid on;
plot(ts, rad2deg(Xs(2,:)), ts, rad2deg(Xs(3,:)), 'LineWidth', 1.2);
plot(tp, rad2deg(Xp(2,:)), 'k:', tp, rad2deg(Xp(3,:)), 'k:');
xline(T, 'k--', 'LQR haut', 'HandleVisibility', 'off');
ylabel('\theta [deg]'); legend('\theta_1', '\theta_2', 'planifié', 'Location', 'best');
title(sprintf('Swing-up du double pendule : TVLQR (Fmax = %g N, perturbation %g°)', Fmax, dth0));

subplot(3,1,2); hold on; grid on;
plot(ts, Xs(1,:), 'LineWidth', 1.2); plot(tp, Xp(1,:), 'k:');
xline(T, 'k--', 'HandleVisibility', 'off'); ylabel('x [m]');

subplot(3,1,3); hold on; grid on;
plot(ts, Fs, 'LineWidth', 1.2); plot(tp, Fplan(tp), 'k:');
xline(T, 'k--', 'HandleVisibility', 'off');
ylabel('F [N]'); xlabel('t [s]'); legend('appliquée', 'planifiée', 'Location', 'best');

%% 8) Animation (GIF)
faire_gif = true;
if faire_gif
    gif = '11_anim_swingup_double.gif';
    fps = 25;  pas = round(1/(fps*dt));
    L   = p(3+N+1 : 3+2*N);
    xmin = min(Xs(1,:)) - 0.4;  xmax = max(Xs(1,:)) + 0.4;
    f  = figure('Color', 'w', 'Position', [100 100 760 480]);
    ax = axes(f);
    premier = true;
    for k = 1:pas:n
        x = Xs(1,k);  th = Xs(2:nq,k);
        px = [x, x + L(1)*sin(th(1)), x + L(1)*sin(th(1)) + L(2)*sin(th(2))];
        py = [0,     L(1)*cos(th(1)),     L(1)*cos(th(1)) + L(2)*cos(th(2))];

        cla(ax); hold(ax, 'on');
        plot(ax, [xmin xmax], [0 0], 'k-', 'LineWidth', 1);
        rectangle(ax, 'Position', [x-0.12, -0.05, 0.24, 0.10], ...
                  'FaceColor', [0.3 0.5 0.8], 'Curvature', 0.2);
        coul = [0.85 0.33 0.1];  if Ms(k) == 2, coul = [0.2 0.6 0.2]; end
        plot(ax, px, py, 'o-', 'LineWidth', 4, 'Color', coul, ...
             'MarkerFaceColor', 'k', 'MarkerSize', 5);
        axis(ax, 'equal');  xlim(ax, [xmin xmax]);  ylim(ax, [-sum(L)-0.1, sum(L)+0.1]);
        txt = 'suivi TVLQR';  if Ms(k) == 2, txt = 'LQR haut'; end
        title(ax, sprintf('t = %.2f s   |   %s   |   F = %+.1f N', ts(k), txt, Fs(k)));
        set(ax, 'YTick', []);  xlabel(ax, 'x [m]');  box(ax, 'on');  drawnow;

        [im, map] = rgb2ind(frame2im(getframe(f)), 256);
        if premier
            imwrite(im, map, gif, 'gif', 'LoopCount', Inf, 'DelayTime', 1/fps);
            premier = false;
        else
            imwrite(im, map, gif, 'gif', 'WriteMode', 'append', 'DelayTime', 1/fps);
        end
    end
    fprintf('GIF enregistré : %s\n', gif);
end

%% ---------- Fonctions locales (tout passe par les arguments) ----------
function dX = dyn(X, F, p, N)
% Modèle non linéaire complet
    nq = N + 1;
    q  = X(1:nq);  dq = X(nq+1:end);
    Q  = [F - p(2)*dq(1); zeros(N,1)];
    dX = [dq; nlink_M(q, p) \ (Q - nlink_h(q, dq, p))];
end

function [A, B] = linearise(X, F, p, N)
% Linéarisation numérique (différences finies centrées) autour de (X, F)
    ns = numel(X);  e = 1e-6;
    A  = zeros(ns);
    for i = 1:ns
        d = zeros(ns,1);  d(i) = e;
        A(:,i) = (dyn(X+d, F, p, N) - dyn(X-d, F, p, N)) / (2*e);
    end
    B = (dyn(X, F+e, p, N) - dyn(X, F-e, p, N)) / (2*e);
end

function ds = riccati(tt, s, tg, dtg, Ag, Bg, Q, R, ns)
% Équation de Riccati différentielle :  -dS/dt = A'S + SA - S B R^-1 B' S + Q
    k = min(max(floor(tt/dtg) + 1, 1), numel(tg) - 1);
    w = (tt - tg(k)) / dtg;
    A = (1-w)*Ag(:,:,k) + w*Ag(:,:,k+1);
    B = (1-w)*Bg(:,k)   + w*Bg(:,k+1);
    Sm = reshape(s, ns, ns);
    dS = -(A.'*Sm + Sm*A - Sm*(B*B.')*Sm/R + Q);
    dS = (dS + dS.')/2;                      % garde S symétrique
    ds = dS(:);
end
