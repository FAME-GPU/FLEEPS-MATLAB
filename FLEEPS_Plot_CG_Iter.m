%% plot CG iter

clear
clc
close all

%% plot FCC
load DATA/LS_info_FCC.mat
n = size(LS_info_SVD, 2);
grid_nums = (1 : n)' * 10;
matrix_size = 3 * grid_nums.^3;
figure
% plot(grid_nums, LS_info_WSVD(1, :), 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'red'), hold on
% plot(grid_nums, LS_info_SVD(1, :), 'bs-', 'LineWidth', 2, 'MarkerFaceColor', 'blue'), hold on
semilogx(matrix_size, LS_info_WSVD(1, :), 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'red'), hold on
semilogx(matrix_size, LS_info_SVD(1, :), 'bd-', 'LineWidth', 2, 'MarkerFaceColor', 'blue'), hold on
xlim([matrix_size(1)/1.2 1.2*matrix_size(end)])
ylim([10 90])
xlabel('matrix size $3n$','FontSize',16,'Interpreter','latex'); 
ylabel('#iter(CG)','FontSize',16);
legend('WSVD-CG', 'SVD-CG', 'Location','east', 'FontSize', 14);
grid on
box off
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on
% h2 = axes('Position', [0.65 0.3 0.25 0.25]);
% plot(matrix_size, CD_FCC, '+', 'Color', 'm', 'LineWidth', 2, 'MarkerFaceColor', 'm'), hold on
% axis(h2, [10 170 min(CD_FCC) - 10 max(CD_FCC) + 10])
% saveas(gca, 'CG_iter_FCC', 'epsc');

%% plot SC
load DATA/LS_info_SC.mat
n = size(LS_info_SVD, 2);
grid_nums = (1 : n)' * 10;
matrix_size = 3 * grid_nums.^3;
figure
% plot(grid_nums, LS_info_WSVD(1, :), 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'red'), hold on
% plot(grid_nums, LS_info_SVD(1, :), 'bs-', 'LineWidth', 2, 'MarkerFaceColor', 'blue'), hold on
semilogx(matrix_size, LS_info_WSVD(1, :), 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'red'), hold on
semilogx(matrix_size, LS_info_SVD(1, :), 'bd-', 'LineWidth', 2, 'MarkerFaceColor', 'blue'), hold on
xlim([matrix_size(1)/1.2 1.2*matrix_size(end)])
ylim([10 90])
xlabel('matrix size $3n$','FontSize',16,'Interpreter','latex'); 
ylabel('#iter(CG)','FontSize',16);
legend('WSVD-CG', 'SVD-CG', 'Location','east', 'FontSize', 14);
grid on
box off
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on
% saveas(gca, 'CG_iter_SC', 'epsc');

%% plot HEX
load DATA/LS_info_HEX.mat
n = size(LS_info_SVD, 2);
grid_nums = (1 : n)' * 10;
matrix_size = 3 * grid_nums.^3;
figure
% plot(grid_nums, LS_info_WSVD(1, :), 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'red'), hold on
% plot(grid_nums, LS_info_SVD(1, :), 'bs-', 'LineWidth', 2, 'MarkerFaceColor', 'blue'), hold on
semilogx(matrix_size, LS_info_WSVD(1, :), 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'red'), hold on
semilogx(matrix_size, LS_info_SVD(1, :), 'bd-', 'LineWidth', 2, 'MarkerFaceColor', 'blue'), hold on
xlim([matrix_size(1)/1.2 1.2*matrix_size(end)])
ylim([10 90])
xlabel('matrix size $3n$','FontSize',16,'Interpreter','latex'); 
ylabel('#iter(CG)','FontSize',16);
legend('WSVD-CG', 'SVD-CG', 'Location','east', 'FontSize', 14);
grid on
box off
ax2 = axes('Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on
% saveas(gca, 'CG_iter_HEX', 'epsc');




