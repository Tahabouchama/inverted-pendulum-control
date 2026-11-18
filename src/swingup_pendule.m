%% swingup_pendule.m — Swing-up par énergie + bascule sur le LQR
% Le pendule part de la position BASSE et doit arriver en position haute.
% Deux commandes qui se relaient :
%   - loin du haut  : pompage d'énergie (loi non linéaire)
%   - près du haut  : LQR (stabilisation)
%
% Prérequis : avoir lancé derive_eom.m (crée cartpole_dyn.m et cartpole_lin.m).
% Requiert : Control System Toolbox.

clear; clc; close all;

%% 1) Paramètres (identiques aux autres scripts)
M = 0.5;  m = 0.2;  l = 0.3;
J = m*(2*l)^2/12;
b = 0.1;  g = 9.81;
p = [M m l J b g];

%% 2) Gain LQR pour la phase finale
[A, B] = cartpole_lin(p);
K = lqr(A, B, diag([10 100 1 1]), 1);

%% 3) Énergie du pendule
Jt   = J + m*l^2;                              % inertie autour du pivot
Eref = m*g*l;                                  % énergie visée : en haut, immobile

%% 4) Réglages du swing-up
kE    = 6;                % gain de pompage (plus grand = montée plus rapide)
kx    = 1.0;              % recentrage du chariot pendant le pompage
kdx   = 2.0;
thSw  = deg2rad(25);      % bascule vers le LQR sous cet angle...
dthSw = 2.0;              % ...et sous cette vitesse angulaire [rad/s]
Fmax  = 10;               % saturation moteur [N]

T = 10;  dt = 1e-3;                            % durée [s] et pas de temps

%% 5) Simulation (Runge-Kutta 4, pas fixe : la commande commute, donc ode45
%    serait ralenti par les discontinuités)
t  = 0:dt:T;
n  = numel(t);
Xs = zeros(n, 4);
Fs = zeros(n, 1);
Ms = zeros(n, 1);

Xs(1,:) = [0, pi, 0, 0];                       % départ : pendule EN BAS, immobile

for k = 1:n
    % --- état courant, angle ramené dans [-pi, pi] ---
    x = Xs(k,1);  dx = Xs(k,3);  dth = Xs(k,4);
    th = mod(Xs(k,2) + pi, 2*pi) - pi;
    E  = 0.5*Jt*dth^2 + m*g*l*cos(th);

    % --- choix de la commande ---
    if abs(th) < thSw && abs(dth) < dthSw
        Ms(k) = 2;                                     % LQR
        F = -K * [x; th; dx; dth];
    else
        Ms(k) = 1;                                     % pompage d'énergie
        s = sign(dth*cos(th));
        if s == 0, s = 1; end
        a = kE*(E - Eref)*s - kx*x - kdx*dx;           % accélération visée du chariot
        % force exacte qui produit cette accélération (modèle complet)
        ddth = (m*g*l*sin(th) - m*l*cos(th)*a) / Jt;
        F = (M+m)*a + m*l*cos(th)*ddth - m*l*sin(th)*dth^2 + b*dx;
    end
    Fs(k) = max(min(F, Fmax), -Fmax);

    % --- intégration RK4, commande maintenue constante sur le pas ---
    if k < n
        X  = Xs(k,:).';
        f1 = cartpole_dyn(X,           Fs(k), p);
        f2 = cartpole_dyn(X + dt/2*f1, Fs(k), p);
        f3 = cartpole_dyn(X + dt/2*f2, Fs(k), p);
        f4 = cartpole_dyn(X + dt*f3,   Fs(k), p);
        Xs(k+1,:) = (X + dt/6*(f1 + 2*f2 + 2*f3 + f4)).';
    end
end

th  = mod(Xs(:,2) + pi, 2*pi) - pi;
E   = 0.5*Jt*Xs(:,4).^2 + m*g*l*cos(th);
kSw = find(Ms == 2, 1);                        % instant de bascule

if isempty(kSw)
    fprintf('Le LQR n''a jamais pris la main : augmenter T ou kE.\n');
else
    fprintf('Bascule vers le LQR à t = %.2f s\n', t(kSw));
end
fprintf('Angle final : %+.2f°   |   Force max : %.1f N\n', ...
        rad2deg(th(end)), max(abs(Fs)));

%% 6) Tracés
figure('Color','w','Name','Swing-up','Position',[100 100 820 700]);

subplot(4,1,1); plot(t, rad2deg(th)); grid on; ylabel('\theta [deg]');
title('Swing-up par énergie puis stabilisation LQR');
subplot(4,1,2); plot(t, Xs(:,1)); grid on; ylabel('x [m]');
subplot(4,1,3);
plot(t, E, t([1 end]), [Eref Eref], '--'); grid on; ylabel('E [J]');
legend('énergie du pendule', 'E_{ref} (position haute)', 'Location','southeast');
subplot(4,1,4); plot(t, Fs); grid on; ylabel('F [N]'); xlabel('t [s]');

if ~isempty(kSw)                               % marque la bascule sur chaque panneau
    for i = 1:4
        subplot(4,1,i); hold on; xline(t(kSw), 'k--', 'LQR', 'HandleVisibility', 'off');
    end
end

%% 7) Animation (GIF)
faire_gif = false;
if faire_gif
    gif = '10_anim_swingup.gif';
    fps = 25;  pas = round(1/(fps*dt));        % on ne garde qu'une image sur "pas"
    xmin = min(Xs(:,1)) - 0.4;  xmax = max(Xs(:,1)) + 0.4;
    L = 2*l;

    f = figure('Color','w','Position',[100 100 760 420]);
    premier = true;
    for k = 1:pas:n
        x = Xs(k,1);  a = th(k);
        clf; hold on;
        plot([xmin xmax], [0 0], 'k-', 'LineWidth', 1);
        rectangle('Position', [x-0.12, -0.05, 0.24, 0.10], ...
                  'FaceColor', [0.3 0.5 0.8], 'Curvature', 0.2);
        coul = [0.85 0.33 0.1];  if Ms(k) == 2, coul = [0.2 0.6 0.2]; end
        plot([x, x + L*sin(a)], [0, L*cos(a)], 'o-', 'LineWidth', 4, ...
             'Color', coul, 'MarkerFaceColor', 'k', 'MarkerSize', 5);
        axis equal; xlim([xmin xmax]); ylim([-L-0.2, L+0.2]);
        txt = 'pompage'; if Ms(k) == 2, txt = 'LQR'; end
        title(sprintf('t = %.2f s   |   %s   |   E = %+.2f J   |   F = %+.1f N', ...
                      t(k), txt, E(k), Fs(k)));
        set(gca,'YTick',[]); xlabel('x [m]'); box on; drawnow;

        [im, map] = rgb2ind(frame2im(getframe(f)), 256);
        if premier
            imwrite(im, map, gif, 'gif', 'LoopCount', Inf, 'DelayTime', 1/fps);
            premier = false;
        else
            imwrite(im, map, gif, 'gif', 'WriteMode','append', 'DelayTime', 1/fps);
        end
    end
    fprintf('GIF enregistré : %s\n', gif);
end
