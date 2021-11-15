%% LOADING THE VARIABLES FOR THE TEST

clear 
clc

load('test_functions.mat')



%% RUN THE STEEPEST DESCENT

disp('**** STEEPEST DESCENT: START *****')

[xk, fk, gradfk_norm, k, xseq] = steepest_descent(x0, f, ...
    gradf, alpha, kmax, tolgrad);
disp('**** STEEPEST DESCENT: FINISHED *****')
disp('**** STEEPEST DESCENT: RESULTS *****')
disp('************************************')
disp(['xk: ', mat2str(xk), ' (actual minimum: [0; 0]);'])
disp(['f(xk): ', num2str(fk), ' (actual min. value: 0);'])
disp(['N. of Iterations: ', num2str(k),'/',num2str(kmax), ';'])
disp('************************************')


%% PLOTS

% Creation of the meshgrid for the contour-plot
[X, Y] = meshgrid(linspace(-6, 6, 500), linspace(-6, 6, 500));
% Computation of the values of f for each point of the mesh
Z = X.^2 + Y.^2;

% Plots

% Simple Plot
fig1 = figure();
% Contour plot with curve levels chosen authomatically
[C1, ~] = contour(X, Y, Z);
hold on
% plot of the points in xseq
plot(xseq(1, :), xseq(2, :), '--*')
hold off

% More interesting Plot
fig2 = figure();
% Contour plot with curve levels for each point in xseq
% ATTENTION: actually, the mesh [X, Y] is to coarse for plotting the last
% level curves corresponding to the last point in xseq (check it zooming
% the image).
[C2, ~] = contour(X, Y, Z, f(xseq));
hold on
% plot of the points in xseq
plot(xseq(1, :), xseq(2, :), '--*')
hold off






