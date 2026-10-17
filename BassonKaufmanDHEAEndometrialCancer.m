% =========================================================================
% BassonKaufmanDHEAEndometrialCancer.m
% Non-Linear Growth Dynamics of Endometrial Adenocarcinoma Under Exogenous
% Influx of DHEA - Biomathematics 374, Stellenbosch University (2026)
% Authors: Helen Basson & Kaylee Kaufman
%
% What it does: solves the 3-ODE model below with ode45, computes the
% critical oestradiol level (E_crit) and critical DHEA dose (S_crit), and
% produces all result figures (saved as Fig1.png ... Fig7.png).
%
% How to run: open in MATLAB (R2019b or later, no toolboxes needed),
% set the current folder to this file's folder, and press Run.
% =========================================================================
%
% Model:
% dD_S/dt = q*S - deltaD*D_S
% dE/dt   = Phi*(Db + D_S) - deltaE*E
% dN/dt   = r*H(E)*N*(1 - N/L) - dT*N,   H(E) = E^n/(Km^n + E^n)

clear; clc; close all;

% ---- 1. PARAMETERS (Table 1 of the report) ----
Db     = 6.31;    % baseline tissue DHEA (nM) - postmenopausal serum DHEA 1.82 ng/mL, Stanczyk et al. 2009
q      = 4.01;    % dose-to-tissue DHEA coefficient (nM/mg) - Stanczyk et al. 2009
deltaD = 16.64;   % DHEA clearance (1/day) = ln2/(60 min half-life), Longcope 1998
deltaE = 3.54;    % oestradiol loss (1/day) = ln2/(4.7 h half-life), Hata et al. 1987
Phi    = 0.0354;  % DHEA->oestradiol conversion (1/day) = 0.01*deltaE, assumed 1% net conversion
r      = 0.704;   % max proliferation (1/day) = ln2/1.25 day doubling + dT, Holinka 1986; Nishida 1985
dT     = 0.15;    % tumour cell death rate (1/day) - assumed
Km     = 1.0;     % half-saturation oestradiol, K (nM) - Holinka et al. 1986
n      = 2;       % Hill coefficient (dimensionless) - idealised ERa dimerisation, Notides et al. 1981
L      = 1e7;     % carrying capacity (cells) - chosen scaling constant
N0     = 5e5;     % initial tumour size (cells) - 5% of L

%% 2. MODEL AND THRESHOLDS
H     = @(E) E.^n ./ (Km^n + E.^n);
Estar = @(S) (Phi/deltaE) * (Db + q*S/deltaD);      % steady-state oestradiol for dose S
model = @(t,y,S) [ q*S - deltaD*y(1); ...           % y = [D_S; E; N]
                   Phi*(Db + y(1)) - deltaE*y(2); ...
                   r*H(y(2))*y(3)*(1 - y(3)/L) - dT*y(3) ];
y0    = [0; Phi*Db/deltaE; N0];                     % no supplement yet

Ecrit = Km * (dT/(r - dT))^(1/n);
Scrit = (deltaD/q) * (deltaE*Ecrit/Phi - Db);
Ebase = Phi*Db/deltaE;
fprintf('E_crit = %.4f nM,  S_crit = %.2f mg/day,  baseline E = %.4f nM\n', Ecrit, Scrit, Ebase);

%FIGURE 1 - HORMONE DYNAMICS (Section 4.1)
doses = [50, Scrit, 2*Scrit];
names = {'S = 50 mg/day', 'S = S_{crit}', 'S = 2 S_{crit}'};
col   = {'b', [0.9 0.5 0.1], 'r'};

figure(1)
for k = 1:3
    [t, y] = ode45(@(t,y) model(t,y,doses(k)), [0 3], y0);
    subplot(1,2,1); hold on; plot(t, Db + y(:,1), 'Color', col{k}, 'LineWidth', 2);
    subplot(1,2,2); hold on; plot(t, y(:,2),      'Color', col{k}, 'LineWidth', 2);
end
subplot(1,2,1); xlabel('Time (days)'); ylabel('Total DHEA (nM)'); title('(a) DHEA');
legend(names, 'Location', 'southeast'); grid on;
subplot(1,2,2); plot([0 3], [Ecrit Ecrit], 'k--'); plot([0 3], [Ebase Ebase], 'k:');
xlabel('Time (days)'); ylabel('Oestradiol (nM)'); title('(b) Oestradiol'); grid on;
text(0.05, Ecrit, 'E_{crit}', 'VerticalAlignment', 'bottom');
text(0.05, Ebase, 'baseline', 'VerticalAlignment', 'bottom');
sgtitle('Hormone levels after a DHEA dose is started');

%FIGURE 2 - TUMOUR DYNAMICS (Section 4.2)
mult  = [0.8 1.0 1.2 1.4];          % doses as multiples of S_crit
style = {'b-', 'c--', 'k:', 'r-'};

figure(2); hold on
for k = 1:4
    S = mult(k)*Scrit;
    [t, y] = ode45(@(t,y) model(t,y,S), [0 160], y0);
    h(k) = plot(t, y(:,3), style{k}, 'LineWidth', 2);
    lab{k} = sprintf('S = %.1f mg/day', S);
    fprintf('S = %.1f mg/day: N(160 days) = %.2e cells,  predicted N* = %.2e\n', ...
            S, y(end,3), L*max(0, 1 - dT/(r*H(Estar(S)))));
end
plot([0 160], [L L], 'k--'); text(2, 0.96*L, 'Carrying capacity (L)');
xlabel('Time (days)'); ylabel('Tumour cells, N(t)');
title('Tumour growth at doses around S_{crit}');
legend(h, lab, 'Location', 'northeast'); grid on; axis([0 160 0 1.1*L]);

% FIGURE 3 GRAPHICAL CHECK OF E_crit AND S_crit (Section 4.2)
E = linspace(0, 1.5, 1000);   S = linspace(0, 400, 1000);

figure(3)
subplot(1,2,1); hold on
plot(E, r*H(E), 'b-', 'LineWidth', 2.5);
plot([0 1.5], [dT dT], 'r--', 'LineWidth', 2);
plot([Ecrit Ecrit], [0 dT], 'k:', 'HandleVisibility', 'off'); plot(Ecrit, dT, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
xlabel('Oestradiol, E (nM)'); ylabel('Rate (1/day)'); title('(a) Against oestradiol');
legend('Proliferation r H(E)', sprintf('Death d_T = %.2f', dT), sprintf('E_{crit} = %.3f nM', Ecrit), 'Location', 'northwest');
axis([0 1.5 0 0.8]); grid on;

subplot(1,2,2); hold on
plot(S, r*H(Estar(S)), 'b-', 'LineWidth', 2.5);
plot([0 400], [dT dT], 'r--', 'LineWidth', 2);
plot([Scrit Scrit], [0 dT], 'k:', 'HandleVisibility', 'off'); plot(Scrit, dT, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
xlabel('DHEA dose, S (mg/day)'); ylabel('Rate (1/day)'); title('(b) Against DHEA dose');
legend('Proliferation r H(E*(S))', sprintf('Death d_T = %.2f', dT), sprintf('S_{crit} = %.1f mg/day', Scrit), 'Location', 'northwest');
axis([0 400 0 0.8]); grid on;
sgtitle('Where proliferation balances cell death');

%% 6. FIGURE 4 - BIFURCATION DIAGRAM (Section 4.3)
Nraw  = 1 - dT ./ (r*H(Estar(S)));   % N*/L, negative below S_crit
below = S <= Scrit;

figure(4); hold on
patch([25 100 100 25], [-0.2 -0.2 1 1], [0.9 0.9 0.9], 'EdgeColor', 'none');   % typical supplement doses
h1 = plot(S(below),  zeros(1,sum(below)),  'b-',  'LineWidth', 2.5);
h2 = plot(S(~below), zeros(1,sum(~below)), 'r--', 'LineWidth', 2.5);
plot(S(~below), Nraw(~below), 'b-',  'LineWidth', 2.5);
plot(S(below),  Nraw(below),  'r--', 'LineWidth', 2.5);
plot([Scrit Scrit], [-0.2 1], 'k:'); text(Scrit + 5, -0.17, sprintf('S_{crit} = %.0f', Scrit));
xlabel('DHEA dose, S (mg/day)'); ylabel('Equilibrium tumour size, N^*/L');
title('Bifurcation diagram');
legend([h1 h2], 'Stable', 'Unstable', 'Location', 'northwest');
axis([0 400 -0.2 1]); grid on; box on;

%% 7. FIGURES 5 AND 6 - SENSITIVITY OF S_crit (Section 4.4)
calcScrit = @(Phi_, K_, d_) (deltaD/q) * ((deltaE ./ Phi_) .* K_ .* (d_ ./ (r - d_)).^(1/n) - Db);

% One parameter at a time, 50% to 150% of baseline
scale = linspace(0.5, 1.5, 300);
figure(5); hold on
plot(scale, calcScrit(scale*Phi, Km, dT), 'LineWidth', 2.5);
plot(scale, calcScrit(Phi, scale*Km, dT), 'LineWidth', 2.5);
plot(scale, calcScrit(Phi, Km, scale*dT), 'LineWidth', 2.5);
plot(1, Scrit, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
xlabel('Parameter relative to baseline'); ylabel('S_{crit} (mg/day)');
title('Sensitivity of S_{crit} to \Phi, K and d');
legend('\Phi', 'K', 'd', 'Baseline', 'Location', 'best'); axis([0.5 1.5 50 450]); grid on; box on;

% K and Phi together
[Kg, Pg] = meshgrid(linspace(0.5*Km, 1.5*Km, 200), linspace(0.5*Phi, 1.5*Phi, 200));
Sg = calcScrit(Pg, Kg, dT);
figure(6)
contourf(Kg, Pg, Sg, 20); colorbar; hold on
[C, hc] = contour(Kg, Pg, Sg, [50 100 150 190 250 300 400], 'k', 'LineWidth', 1);
clabel(C, hc);
plot(Km, Phi, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 9); text(Km + 0.03, Phi, 'Baseline');
xlabel('K (nM)'); ylabel('\Phi (day^{-1})');
title('S_{crit} when K and \Phi change together'); box on;

%% 8. SECTION 3 CHECKS WITH ODE SOLVERS
opts  = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
atSS  = @(S, N) [q*S/deltaD; Estar(S); N];          % hormones already at steady state

% 3.1 Equilibria: run the full model for a long time and compare with the formulas
fprintf('\n3.1 Equilibria (ODE vs formula)\n');
for S = [50, 0.8*Scrit, 1.2*Scrit, 2*Scrit]
    [~, y] = ode45(@(t,y) model(t,y,S), [0 1000], y0, opts);
    Nform = L*max(0, 1 - dT/(r*H(Estar(S))));
    fprintf('S = %6.1f:  D_S %.3f vs %.3f   E %.4f vs %.4f   N %.3e vs %.3e\n', ...
            S, y(end,1), q*S/deltaD, y(end,2), Estar(S), y(end,3), Nform);
end

% 3.2 Stability: eigenvalues of a numerical Jacobian vs -deltaD, -deltaE and lambda_3
fprintf('\n3.2 Eigenvalues (numerical vs formula)\n');
for S = [0.8*Scrit, 1.2*Scrit]
    Nstar = L*max(0, 1 - dT/(r*H(Estar(S))));
    for Neq = unique([0, Nstar])
        yeq = atSS(S, Neq);
        J = zeros(3);
        for j = 1:3
            dy = zeros(3,1); dy(j) = 1e-6*max(1, abs(yeq(j)));
            J(:,j) = (model(0, yeq+dy, S) - model(0, yeq-dy, S)) / (2*dy(j));
        end
        lam3 = r*H(Estar(S))*(1 - 2*Neq/L) - dT;
        fprintf('S = %6.1f, N* = %.2e:  numerical [%s]   formula [%.4f %.4f %.4f]\n', ...
                S, Neq, sprintf('%.4f ', sort(eig(J))), sort([-deltaD -deltaE lam3]));
    end
end

% 3.3 and 3.4: find S_crit where a small tumour neither grows nor shrinks
Scrit_ode = fzero(@(S) growthRate(model, S, atSS(S, 1e3), 50, opts), [0.5 1.5]*Scrit);
[~, y]    = ode45(@(t,y) model(t,y,Scrit_ode), [0 10], y0, opts);
fprintf('\n3.3 E_crit: ODE %.4f vs formula %.4f nM\n', y(end,2), Ecrit);
fprintf('3.4 S_crit: ODE %.2f vs formula %.2f mg/day\n', Scrit_ode, Scrit);

%% 9. SAVE ALL FIGURES (optional)
for k = 1:6, saveas(figure(k), sprintf('Fig%d.png', k)); end

%% LOCAL FUNCTION
function g = growthRate(model, S, yStart, T, opts)
    % average growth rate of a small tumour over T days
    [~, y] = ode45(@(t,y) model(t,y,S), [0 T], yStart, opts);
    g = log(y(end,3)/yStart(3)) / T;
end
