%% derive_eom_nlink.m — Pendule inversé à N barres en série sur chariot
% Généralisation de derive_eom.m. Modèle de Lagrange sous forme matricielle :
%       M(q)*qdd + h(q,qd) = Q
% puis linéarisation autour de la position tout-debout.
%
% Conventions :
%   q   = [x; th1; ...; thN]   x : chariot [m], thi : angle ABSOLU de la barre i
%                              depuis la verticale haute [rad] (0 = debout)
%   Barre i : tige homogène, masse mi, longueur Li, CdG au milieu, J = mi*Li^2/12
%   Etat  X = [q; qd]  (2N+2 états),  commande F [N] sur le chariot
%
% Astuce de validation : avec N = 1, A et B doivent coïncider avec derive_eom.m
% (en prenant l = L/2 et J = m*L^2/12).
%
% Requiert : Symbolic Math Toolbox (Control System Toolbox pour la suite)

clear; clc;
N = 2;                                   % nombre de barres (tester d'abord 1, puis 2)

%% 1) Variables symboliques
syms x dx F real
syms M0 b g positive
th  = sym('th',  [N 1], 'real');
dth = sym('dth', [N 1], 'real');
mL  = sym('m',   [N 1], 'positive');
Ls  = sym('L',   [N 1], 'positive');

q  = [x; th];
dq = [dx; dth];
p_sym = [M0; b; g; mL; Ls];              % vecteur de paramètres

%% 2) Énergies, barre par barre
T = 1/2*M0*dx^2;
V = sym(0);
px = x; py = sym(0);                     % articulation courante (départ : pivot sur le chariot)
for i = 1:N
    cx = px + Ls(i)/2*sin(th(i));        % centre de gravité de la barre i
    cy = py + Ls(i)/2*cos(th(i));
    vc = jacobian([cx; cy], q) * dq;     % vitesse du CdG
    Ji = mL(i)*Ls(i)^2/12;
    T  = T + 1/2*mL(i)*(vc.'*vc) + 1/2*Ji*dth(i)^2;
    V  = V + mL(i)*g*cy;
    px = px + Ls(i)*sin(th(i));          % articulation suivante
    py = py + Ls(i)*cos(th(i));
end

%% 3) Forme matricielle  M(q)*qdd + h(q,qd) = Q
Mq = simplify(hessian(T, dq));                         % matrice de masse (symétrique, définie positive)
h  = jacobian(Mq*dq, q)*dq - jacobian(T - V, q).';     % Coriolis/centrifuge + gravité

Bu = [1; zeros(N,1)];                                  % la force n'agit que sur x
D  = diag([b; zeros(N,1)]);                            % frottement visqueux du chariot
% Q = Bu*F - D*dq

%% 4) Linéarisation autour de q = 0, qd = 0
%   M0m*qdd + K*q + D*qd = Bu*F   (les termes en qd^2 disparaissent)
M0m = simplify(subs(Mq, q, zeros(N+1,1)));
K   = simplify(subs(jacobian(h, q), [q; dq], zeros(2*N+2,1)));

%% 5) Génération des fonctions numériques
matlabFunction(Mq,     'File', 'nlink_M',    'Vars', {q, p_sym});
matlabFunction(h,      'File', 'nlink_h',    'Vars', {q, dq, p_sym});
matlabFunction(M0m, K, 'File', 'nlink_linM', 'Vars', {p_sym});
fprintf('Fichiers générés : nlink_M.m, nlink_h.m, nlink_linM.m\n');

%% 6) Paramètres numériques
M0v = 1.0;  bv = 0.1;  gv = 9.81;
mv  = 0.10*ones(N,1);                    % kg par barre
Lv  = 0.30*ones(N,1);                    % m par barre
p   = [M0v; bv; gv; mv; Lv];

%% 7) Modèle d'état linéaire  dX = A*X + B*F
[M0n, Kn] = nlink_linM(p);
Dn  = diag([bv; zeros(N,1)]);
nq  = N+1;
A = [zeros(nq), eye(nq);
     -M0n\Kn,   -M0n\Dn];
B = [zeros(nq,1);
     M0n\Bu];

fprintf('\nPôles en boucle ouverte (N = %d) :\n', N); disp(eig(A))
fprintf('Nombre de pôles instables : %d (attendu : %d)\n', sum(real(eig(A)) > 1e-9), N);

% Commandabilité (sans toolbox) — attention au conditionnement pour N grand
Co = B; for k = 1:2*nq-1, Co = [Co, A*Co(:,end)]; end %#ok<AGROW>
fprintf('Rang de commandabilité : %d / %d   (cond = %.1e)\n', rank(Co), 2*nq, cond(Co));

save('nlink_model.mat', 'A', 'B', 'p', 'N');

%% 8) Simulation non linéaire en boucle ouverte
dyn = @(X, Fin, pp) [X(nq+1:end);
                     nlink_M(X(1:nq), pp) \ ( [1; zeros(N,1)]*Fin ...
                                            - diag([pp(2); zeros(N,1)])*X(nq+1:end) ...
                                            - nlink_h(X(1:nq), X(nq+1:end), pp) )];

X0 = [0; 0.01*ones(N,1); zeros(nq,1)];
[t, Xs] = ode45(@(t,X) dyn(X, 0, p), [0 2], X0);

figure('Name', 'Boucle ouverte N barres');
plot(t, rad2deg(Xs(:,2:nq))); grid on;
xlabel('t [s]'); ylabel('\theta_i [deg]');
legend(compose('\\theta_%d', 1:N)); title('Boucle ouverte : effondrement (attendu)');

%% 9) Contrôle de cohérence : conservation de l'énergie avec b = 0
p0 = p; p0(2) = 0;
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
[t0, X0s] = ode45(@(t,X) dyn(X, 0, p0), [0 1], X0, opts);
Efun = matlabFunction(T + V, 'Vars', {q, dq, p_sym});
E = arrayfun(@(k) Efun(X0s(k,1:nq).', X0s(k,nq+1:end).', p0), 1:numel(t0));
fprintf('Dérive d''énergie max (b = 0) : %.2e J\n', max(abs(E - E(1))));
