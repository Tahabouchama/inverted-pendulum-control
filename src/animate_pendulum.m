%% animate_pendulum.m — Animation du pendule inversé (simple ou double)
% Rejoue une simulation image par image et enregistre un GIF dans ../figures/.
% Prérequis : avoir lancé derive_eom.m et derive_eom_nlink.m (N = 2).
% À lancer depuis le dossier src/.

clear; clc; close all;

%% ----- Réglages -----
systeme = 'simple';   % 'simple' ou 'double'
angle0  = 45;         % inclinaison initiale [deg]
Fmax    = 10;         % saturation moteur [N]
T       = 5;          % durée simulée [s]
fps     = 25;         % images par seconde du GIF

%% ----- Modèle, gain et dynamique -----
switch systeme
    case 'simple'
        M = 0.5; m = 0.2; l = 0.3; J = m*(2*l)^2/12; b = 0.1; g = 9.81;
        p = [M m l J b g];
        [A, B] = cartpole_lin(p);
        N = 1;  nq = 2;
        L = 2*l*ones(N,1);                       % longueur de la barre (tige de longueur 2l)
        Q = diag([10 100 1 1]);
        dyn = @(X, F) cartpole_dyn(X, F, p);
        gif = '08_anim_simple.gif';

    case 'double'
        S  = load('nlink_model.mat', 'A', 'B', 'p', 'N');
        A = S.A;  B = S.B;  p = S.p;  N = S.N;  nq = N + 1;
        L = p(3+N+1 : 3+2*N);                    % longueurs des barres
        Q = diag([10, 100*ones(1,N), ones(1,nq)]);
        dyn = @(X, F) [ X(nq+1:end);
                        nlink_M(X(1:nq), p) \ ( [F; zeros(N,1)] ...
                                              - [p(2)*X(nq+1); zeros(N,1)] ...
                                              - nlink_h(X(1:nq), X(nq+1:end), p) ) ];
        gif = '09_anim_double.gif';
end

K   = lqr(A, B, Q, 1);
loi = @(X) max(min(-K*X, Fmax), -Fmax);

%% ----- Simulation, échantillonnée à pas de temps régulier -----
X0 = zeros(2*nq, 1);  X0(2:nq) = deg2rad(angle0);
t  = linspace(0, T, T*fps);
[~, Xs] = ode45(@(t,X) dyn(X, loi(X)), t, X0);

%% ----- Animation -----
Ltot = sum(L);
xmin = min(Xs(:,1)) - 0.4;  xmax = max(Xs(:,1)) + 0.4;

f = figure('Color', 'w', 'Position', [100 100 760 420]);
for k = 1:numel(t)
    x  = Xs(k,1);  th = Xs(k, 2:nq);

    % Positions des articulations, de la base du chariot au sommet
    px = x;  py = 0;
    for i = 1:N
        px(end+1) = px(end) + L(i)*sin(th(i)); %#ok<SAGROW>
        py(end+1) = py(end) + L(i)*cos(th(i)); %#ok<SAGROW>
    end

    clf; hold on;
    plot([xmin xmax], [0 0], 'k-', 'LineWidth', 1);                    % rail
    rectangle('Position', [x-0.12, -0.05, 0.24, 0.10], ...             % chariot
              'FaceColor', [0.3 0.5 0.8], 'Curvature', 0.2);
    plot(px, py, 'o-', 'LineWidth', 4, 'Color', [0.85 0.33 0.1], ...   % barres
         'MarkerFaceColor', 'k', 'MarkerSize', 5);

    axis equal; xlim([xmin xmax]); ylim([-0.2, Ltot + 0.2]);
    title(sprintf('t = %.2f s   |   \\theta_1 = %+.1f°   |   F = %+.1f N', ...
                  t(k), rad2deg(th(1)), loi(Xs(k,:).')));
    set(gca, 'YTick', []); xlabel('x [m]'); box on; drawnow;

    % Ajout de l'image au GIF
    [im, map] = rgb2ind(frame2im(getframe(f)), 256);
    if k == 1
        imwrite(im, map, gif, 'gif', 'LoopCount', Inf, 'DelayTime', 1/fps);
    else
        imwrite(im, map, gif, 'gif', 'WriteMode', 'append', 'DelayTime', 1/fps);
    end
end

fprintf('GIF enregistré : %s\n', gif);
