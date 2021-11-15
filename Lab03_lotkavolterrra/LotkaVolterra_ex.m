%% LOADING THE VARIABLES FOR THE TEST

clear 
close all
clc

load('PoliTO_preypreds_guided.mat')


%% SOLUTION COMPUTED WITH LOTKA-VOLTERRA MODEL

[xn,yn,x_eq,y_eq,tn] = ...
    Lotka_Volterra_E_Espl(a,b,c,d, x0,y0,t0,T,N);

%% VISUALIZATION OF THE SOLUTION COMPUTED WITH LOTKA-VOLTERRA MODEL


fig_sol1 = figure();
plot(xn, 'LineWidth', 2)
hold on
plot(yn, 'LineWidth', 2)
legend('preys', 'predators')
xlabel('time')
ylabel('population size')
fig_sol1.Children(1).FontSize=20;
fig_sol1.Children(2).FontSize=20;
hold off

fig_sol2 = figure();
plot(xn, yn)
hold on
plot(x_eq, y_eq, '*', 'LineWidth', 2)
xlabel('preys')
ylabel('predators')
legend('pop.s periodic evolution', 'stable eq. point')
fig_sol2.Children(1).FontSize=20;
fig_sol2.Children(2).FontSize=20;
hold off


%% START OF THE "FIND-THE-COEFFICIENTS PROBLEM"

% REMEMBER: 
% We are assuming that the values of a,b,c,d,x0,y0 are unkown.
% We are assuming that the only values are the noisy-data inside 
% the table preypred_samples.

prey_pred_samples = preypred_samples{:,2:end};

% VISUALIZATION OF THE NOISY DATA AND ITS DIFFERENT BEHAVIOUR 
% WITH RESPECT TO THE REAL (UNKNOWN) SOLUTION:

fig_comparison_start1 = figure();
plot(prey_pred_samples(:, 1), '--*', 'LineWidth', 2)
hold on
plot(prey_pred_samples(:, 2), '--*', 'LineWidth', 2)
legend('sampled preys', 'sampled predators')
xlabel('years')
ylabel('population size')
fig_comparison_start1.Children(1).FontSize=20;
fig_comparison_start1.Children(2).FontSize=20;
hold off


fig_comparison_start2 = figure();
plot(xn, yn)
hold on
plot(x_eq, y_eq, '*', 'LineWidth', 2)
plot(prey_pred_samples(:, 1), prey_pred_samples(:, 2), ...
    '--*', 'LineWidth', 2)
legend('pop.s periodic evolution', 'stable eq. point', ...
    'population samplings')
fig_comparison_start2.Children(1).FontSize=20;
fig_comparison_start2.Children(2).FontSize=20;
hold off

%% DEFINE THE LOSS FUNCTION AND ITS GRADIENT

% coeffs = [a; b; c; d]  % ATTENTION: IT IS A COLUMN VECTOR!

% DEFINITION OF H
% ATTENTION: prey AND predator ARE COLUMN VECTORS!
H = @(coeffs, prey, predator)...
    coeffs(3) .* log(prey) - coeffs(4) .* prey - coeffs(2) .* predator ...
    + coeffs(1) .* log(predator);

% DEFINITION OF H AS FUNCTION OF coeffs, 
% WITH SAMPLES prey_pred_samples AS PARAMETERS
H_samples = @(coeffs) ...
    H(coeffs, prey_pred_samples(:,1), prey_pred_samples(:,2));

% DEFINITION OF "H-BAR" AS FUNCTION OF coeffs, 
% WITH SAMPLES prey_pred_samples AS PARAMETERS
Hmean_samples = @(coeffs) mean(H_samples(coeffs));

% DEFINITION OF THE LOSS FUNCTION
Loss = @(coeffs) ...
    mean((H_samples(coeffs) - Hmean_samples(coeffs)).^2);

% DEFINITION OF THE MANY PARTIAL DERIVATIVES
H_a_samples = @(coeffs) 0;  % WRITE THE CORRECT EXPRESSION
H_b_samples = @(coeffs) 0;  % WRITE THE CORRECT EXPRESSION
H_c_samples = @(coeffs) 0;  % WRITE THE CORRECT EXPRESSION
H_d_samples = @(coeffs) 0;  % WRITE THE CORRECT EXPRESSION

Hmean_a_samples = @(coeffs) 0;  % WRITE THE CORRECT EXPRESSION
Hmean_b_samples = @(coeffs) 0;  % WRITE THE CORRECT EXPRESSION
Hmean_c_samples = @(coeffs) 0;  % WRITE THE CORRECT EXPRESSION
Hmean_d_samples = @(coeffs) 0;  % WRITE THE CORRECT EXPRESSION

Loss_a = @(coeffs) 0;  % WRITE THE CORRECT EXPRESSION
Loss_b = @(coeffs) 0;  % WRITE THE CORRECT EXPRESSION
Loss_c = @(coeffs) 0;  % WRITE THE CORRECT EXPRESSION
Loss_d = @(coeffs) 0;  % WRITE THE CORRECT EXPRESSION

% DEFINITION OF THE GRADIENT OF THE LOSS FUNCTION
gradLoss = @(coeffs) ...
    [Loss_a(coeffs); Loss_b(coeffs); ...
    Loss_c(coeffs); Loss_d(coeffs)];



%% RUN THE STEEPEST DESCENT

% ... INITIALIZE THE PARAMETERS FOR STEEPEST DESCENT 
% (SEE THE TEXT OF THE EXERCISE)...

coeffs0 = [0; 0; 0; 0]; % WRITE THE CORRECT EXPRESSION
alpha = 0; % WRITE THE CORRECT EXPRESSION
kmax = 0; % WRITE THE CORRECT EXPRESSION
tollgrad = 0; % WRITE THE CORRECT EXPRESSION

disp('**** STEEPEST DESCENT: START *****')

[coeffsk, Lossk, gradLossk_norm, k, xseq] = ...
    steepest_descent(coeffs0, Loss, gradLoss, alpha, kmax, tollgrad);

disp('**** STEEPEST DESCENT: FINISHED *****')
disp('**** STEEPEST DESCENT: RESULTS *****')
disp('************************************')
disp(['pred. a = ', num2str(coeffsk(1)), '(actual a = ', num2str(a), ');'])
% ... DISPLAY THE VALUES ALSO FOR b, c, d ...
disp('****************')
disp(['Loss: ', num2str(Lossk)])
% WE COMPUTE SOME STATS ABOUT THE H VALUES FOR THE SAMPLES,
% USING THE PREDICTED COEFFICIENTS coeffsk
disp(['Hmean_samples: ', num2str(Hmean_samples(coeffsk))])
% THE SMALLER IS THE STANDARD DEVIATION OF H_samples W.R.T. coeffsk,
% THE BETTER THE SOLUTION IS (WE ASSUME), BEACUSE H IS "QUITE-CONSTANT"
% FOR ALL THE SAMPLES.
disp(['Hstd_sampes: ', num2str(std(H_samples(coeffsk)))])
disp('****************')


%% COMPARING SOLUTIONS WITH ACTUAL VALUES (PLOTs)

% Computing populations w.r.t. coeffsk, x0tilde, y0tilde
[xnhat,ynhat,x_eq_hat,y_eq_hat,tn_hat] = ...
    Lotka_Volterra_E_Espl(coeffsk(1),coeffsk(2),...
    coeffsk(3),coeffsk(4),...
    prey_pred_samples(1,1), prey_pred_samples(1,2), ...
    t0,T,N);

% PLOTS W.R.T. TIME (SAMPLES)
fig_comp1b = figure();
plot(tn, xn, 'LineWidth', 2)
hold on
plot(tn, yn, 'LineWidth', 2)
plot(0:size(prey_pred_samples,1)-1, prey_pred_samples(:, 1), ...
    '--*', 'LineWidth', 2)
hold on
plot(0:size(prey_pred_samples,1)-1, prey_pred_samples(:, 2), ...
    '--*', 'LineWidth', 2)
plot(tn, xnhat, 'LineWidth', 2)
plot(tn, ynhat, 'LineWidth', 2)
legend('preys', 'predators', 'prey samples', 'predator samples', ...
    'pred. preys', 'pred. predators')
xlabel('time')
ylabel('population size')
fig_comp1b.Children(1).FontSize=20;
fig_comp1b.Children(2).FontSize=20;
hold off

% PLOTS W.R.T. (prey, predators) PLANE
% Finding the max/min values in the population samples
min_prey = min(prey_pred_samples(:, 1));
max_prey = max(prey_pred_samples(:, 1));
min_pred = min(prey_pred_samples(:, 2));
max_pred = max(prey_pred_samples(:, 2));

% Preparing the Meshgrid for contour plot of function H
% We compute the mesh w.r.t. values greater than the maximums in the
% samplings
[Prey, Pred] = meshgrid(linspace(0, max_prey + 0.4 * max_prey, 100), ...
    linspace(0, max_pred + 0.4 * max_pred, 100));

% Computation of values of H w.r.t. the coefficients computed
% with steepest descent
H_grid = H(coeffsk, Prey, Pred);


fig_comp2 = figure();
plot(xn, yn)
hold on
plot(xnhat, ynhat)
contour(Prey, Pred, H_grid, ...
    [Hmean_samples(coeffsk), Hmean_samples(coeffsk)])
plot(x_eq, y_eq, '*', 'LineWidth', 2)
plot(x_eq_hat, y_eq_hat, '*', 'LineWidth', 2)
xlabel('preys')
ylabel('predators')
legend('pop.s periodic evolution', ...
    'pop.s PREDICTED evolution', ...
    'Hmean_{samples} computed', ...
    'stable eq. point', ...
    'PREDICTED stable eq. point')
fig_comp2.Children(1).FontSize=20;
fig_comp2.Children(2).FontSize=20;
hold off