%% validation_simscape.m — Validation du modèle de Lagrange par Simscape Multibody
% Compare la chute libre du pendule simple calculée par :
%   1) le modèle de Lagrange (cartpole_dyn.m, issu de derive_eom.m)
%   2) le modèle physique Simscape Multibody (pendule_simscape.slx)
% Les deux modèles sont construits indépendamment : s'ils donnent la même
% trajectoire, la dérivation des équations est validée.
%
% Prérequis : derive_eom.m (cartpole_dyn.m) et pendule_simscape.slx
%             dans le dossier courant.
% Requiert  : Simulink, Simscape Multibody.

clear; clc; close all;

mdl = 'pendule_simscape';

%% 1) Paramètres (identiques au modèle Simscape)
M = 0.5;  m = 0.2;  l = 0.3;  J = m*(2*l)^2/12;  b = 0.1;  g = 9.81;
p   = [M m l J b g];
th0 = 0.05;                     % inclinaison initiale [rad], comme dans le Revolute Joint
T   = 3;                        % durée [s]

%% 2) Modèle de Lagrange
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
[tL, XL] = ode45(@(t,X) cartpole_dyn(X, 0, p), [0 T], [0; th0; 0; 0], opts);

%% 3) Modèle Simscape Multibody
load_system(mdl);
set_param(mdl, 'StopTime', num2str(T), 'RelTol', '1e-8', 'AbsTol', '1e-9', 'MaxStep', '1e-3');
out = sim(mdl);

xs  = out.get('x_s');
ths = out.get('th_s');
tS  = xs.Time;
xS  = squeeze(xs.Data);
thS = squeeze(ths.Data);

%% 4) Comparaison sur une grille de temps commune
tc   = linspace(0, T, 1500);
xL   = interp1(tL, XL(:,1), tc);   thL = interp1(tL, XL(:,2), tc);
xSi  = interp1(tS, xS,      tc);   thSi = interp1(tS, thS,    tc);

ex  = abs(xL - xSi);
eth = abs(rad2deg(thL - thSi));
fprintf('Écart max sur x     : %.2e m\n',   max(ex));
fprintf('Écart max sur theta : %.2e deg\n', max(eth));
fprintf('Amplitude de theta  : %.0f deg (le pendule fait de grands balancements)\n', ...
        max(abs(rad2deg(thL))));

%% 5) Tracés
figure('Color', 'w', 'Name', 'Validation Simscape', 'Position', [100 100 820 700]);

subplot(3,1,1); hold on; grid on;
plot(tc, rad2deg(thL), 'LineWidth', 2);
plot(tc, rad2deg(thSi), '--', 'LineWidth', 1.5);
ylabel('\theta [deg]'); legend('Lagrange', 'Simscape Multibody', 'Location', 'best');
title('Validation du modèle : chute libre depuis \theta_0 = 0,05 rad');

subplot(3,1,2); hold on; grid on;
plot(tc, xL, 'LineWidth', 2);
plot(tc, xSi, '--', 'LineWidth', 1.5);
ylabel('x [m]');

subplot(3,1,3);
semilogy(tc, max(eth, 1e-12), tc, max(ex, 1e-12), 'LineWidth', 1.2); grid on;
ylabel('écart'); xlabel('t [s]');
legend('|\Delta\theta| [deg]', '|\Delta x| [m]', 'Location', 'best');
