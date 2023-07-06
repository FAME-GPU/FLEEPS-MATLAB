%% plotConvergenceOrder

clear
clc
close all

load DATA/Eigenvalues_FCC_X.mat
[m, n] = size(eigenvalue);
re = abs(eigenvalue(:, 1 : n - 1) - repmat(eigenvalue(:, end), 1, n - 1)) ./ repmat(eigenvalue(:, end), 1, n - 1);
figure
xx = 1 : n - 1;
yy = log2(re(:, 1 : n - 1));
plot(xx, yy(1, :), 'ro--', 'MarkerSize', 10), hold on
plot(xx, yy(2, :), 'bs--', 'MarkerSize', 10), hold on
plot(xx, yy(3, :), 'kd--', 'MarkerSize', 10), hold on
plot(xx, yy(4, :), 'mp--', 'MarkerSize', 10), hold on 
plot(xx, yy(5, :), 'cd--', 'MarkerSize', 10), hold on
plot(xx, yy(6, :), 'gh--', 'MarkerSize', 10), hold on
plot(xx, yy(7, :), 'Color', [0.8500 0.3250 0.0980], ...
    'Marker', 'x', 'MarkerSize', 10, 'LineStyle', '--'), hold on
idx = 1 : n - 1;
order = zeros(n, 1);
for i = 1 : n
    p = polyfit(xx(idx), yy(i, idx), 1);
    % fprintf(['lambda_' num2str(i) ': %f\n'], -p(1));
    order(i) = -p(1);
end
order_mean = mean(order); 
xx2 = [2 3 4 5];
yy2 = [max(yy(end, :)) + 0.6 0 0 0];
yy2(2 : end) = yy2(1) - order_mean * (1 : length(xx2) - 1);
plot(xx2, yy2, 'k-', 'LineWidth', 2)
text(xx2(2), yy2(2) + 0.5, ['slope = ' num2str(order_mean)], 'FontSize', 16);
set(gca, 'xtick', [1 2 3 4 5 6])
set(gca, 'xticklabel', {'2^{-2}','2^{-3}','2^{-4}','2^{-5}','2^{-6}','2^{-7}'}, 'FontSize', 14)
set(gca, 'ytick', [-16 -14 -12 -10 -8 -6 -4])
set(gca, 'yticklabel', {'2^{-16}','2^{-14}','2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}'}, 'FontSize', 14)
xlabel('$h=1/n_1=1/n_2=1/n_3$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$R_m^h$', 'Interpreter', 'latex', 'FontSize', 20)
% ylabel('$\mathrm{log}_2 (\mathrm{rel\_err})$', 'Interpreter', 'latex', 'FontSize', 20)
legend('$\omega_1$', '$\omega_2$', '$\omega_3$', '$\omega_4$', '$\omega_5$', ...
       '$\omega_6$', '$\omega_7$', ... %'$\omega_8$', '$\omega_9$', '$\omega_{10}$', ... 
       'Location', 'Southwest', 'Interpreter', 'latex', 'FontSize', 14)
axis([1 6 min(min(yy(1 : m, :))) max(max(yy(1 : m, :)))])
grid on

load DATA/Eigenvalues_HEX_A.mat
[m, n] = size(eigenvalue);
re = abs(eigenvalue(:, 1 : n - 1) - repmat(eigenvalue(:, end), 1, n - 1)) ./ repmat(eigenvalue(:, end), 1, n - 1);
figure
xx = 1 : n - 1;
yy = log2(re(:, 1 : n - 1));
plot(xx, yy(1, :), 'ro--', 'MarkerSize', 10), hold on
plot(xx, yy(2, :), 'bs--', 'MarkerSize', 10), hold on
plot(xx, yy(3, :), 'kd--', 'MarkerSize', 10), hold on
plot(xx, yy(4, :), 'mp--', 'MarkerSize', 10), hold on 
plot(xx, yy(5, :), 'cd--', 'MarkerSize', 10), hold on
plot(xx, yy(6, :), 'gh--', 'MarkerSize', 10), hold on
plot(xx, yy(7, :), 'Color', [0.8500 0.3250 0.0980], ...
    'Marker', 'x', 'MarkerSize', 10, 'LineStyle', '--'), hold on
idx = 1 : n - 1;
order = zeros(n, 1);
for i = 1 : n
    p = polyfit(xx(idx), yy(i, idx), 1);
    % fprintf(['lambda_' num2str(i) ': %f\n'], -p(1));
    order(i) = -p(1);
end
order_mean = mean(order); 
xx2 = [2 3 4 5];
yy2 = [max(yy(end, :)) + 0.6 0 0 0];
yy2(2 : end) = yy2(1) - order_mean * (1 : length(xx2) - 1);
plot(xx2, yy2, 'k-', 'LineWidth', 2)
text(xx2(2), yy2(2) + 0.5, ['slope = ' num2str(order_mean)], 'FontSize', 16);
set(gca, 'xtick', [1 2 3 4 5 6])
set(gca, 'xticklabel', {'2^{-2}','2^{-3}','2^{-4}','2^{-5}','2^{-6}','2^{-7}'}, 'FontSize', 14)
set(gca, 'ytick', [-16 -14 -12 -10 -8 -6 -4])
set(gca, 'yticklabel', {'2^{-16}','2^{-14}','2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}'}, 'FontSize', 14)
xlabel('$h=1/n_1=1/n_2=1/n_3$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$R_m^h$', 'Interpreter', 'latex', 'FontSize', 20)
% ylabel('$\mathrm{log}_2 (\mathrm{rel\_err})$', 'Interpreter', 'latex', 'FontSize', 20)
legend('$\omega_1$', '$\omega_2$', '$\omega_3$', '$\omega_4$', '$\omega_5$', ...
       '$\omega_6$', '$\omega_7$', ... %'$\omega_8$', '$\omega_9$', '$\omega_{10}$', ... 
       'Location', 'Southwest', 'Interpreter', 'latex', 'FontSize', 14)
axis([1 6 min(min(yy(1 : m, :))) max(max(yy(1 : m, :)))])
grid on

load DATA/Eigenvalues_SC_M.mat
[m, n] = size(eigenvalue);
re = abs(eigenvalue(:, 1 : n - 1) - repmat(eigenvalue(:, end), 1, n - 1)) ./ repmat(eigenvalue(:, end), 1, n - 1);
figure
xx = 1 : n - 1;
yy = log2(re(:, 1 : n - 1));
plot(xx, yy(1, :), 'ro--', 'MarkerSize', 10), hold on
plot(xx, yy(2, :), 'bs--', 'MarkerSize', 10), hold on
plot(xx, yy(3, :), 'kd--', 'MarkerSize', 10), hold on
plot(xx, yy(4, :), 'mp--', 'MarkerSize', 10), hold on 
plot(xx, yy(5, :), 'cd--', 'MarkerSize', 10), hold on
plot(xx, yy(6, :), 'gh--', 'MarkerSize', 10), hold on
plot(xx, yy(7, :), 'Color', [0.8500 0.3250 0.0980], ...
    'Marker', 'x', 'MarkerSize', 10, 'LineStyle', '--'), hold on
idx = 1 : n - 1;
order = zeros(n, 1);
for i = 1 : n
    p = polyfit(xx(idx), yy(i, idx), 1);
    % fprintf(['lambda_' num2str(i) ': %f\n'], -p(1));
    order(i) = -p(1);
end
order_mean = mean(order); 
xx2 = [2 3 4 5];
yy2 = [max(yy(end, :)) - 0.6 0 0 0];
yy2(2 : end) = yy2(1) - order_mean * (1 : length(xx2) - 1);
plot(xx2, yy2, 'k-', 'LineWidth', 2)
text(xx2(2), yy2(2) + 0.5, ['slope = ' num2str(order_mean)], 'FontSize', 16);
set(gca, 'xtick', [1 2 3 4 5 6])
set(gca, 'xticklabel', {'2^{-2}','2^{-3}','2^{-4}','2^{-5}','2^{-6}','2^{-7}'}, 'FontSize', 14)
set(gca, 'ytick', [-16 -14 -12 -10 -8 -6 -4])
set(gca, 'yticklabel', {'2^{-16}','2^{-14}','2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}'}, 'FontSize', 14)
xlabel('$h=1/n_1=1/n_2=1/n_3$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$R_m^h$', 'Interpreter', 'latex', 'FontSize', 20)
% ylabel('$\mathrm{log}_2 (\mathrm{rel\_err})$', 'Interpreter', 'latex', 'FontSize', 20)
legend('$\omega_1$', '$\omega_2$', '$\omega_3$', '$\omega_4$', '$\omega_5$', ...
       '$\omega_6$', '$\omega_7$', ... %'$\omega_8$', '$\omega_9$', '$\omega_{10}$', ... 
       'Location', 'Southwest', 'Interpreter', 'latex', 'FontSize', 14)
axis([1 6 min(min(yy(1 : m, :))) max(max(yy(1 : m, :)))])
grid on
