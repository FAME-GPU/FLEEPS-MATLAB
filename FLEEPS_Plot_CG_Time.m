%% plot CG time

clear
clc
close all

%% plot FCC
load DATA/LS_info_FCC.mat
n = size(LS_info_SVD, 2);
grid_nums = (1 : n)' * 10;
matrix_size = 3 * grid_nums.^3;
figure
loglog(matrix_size, LS_info_WSVD(2, :), 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'red'), hold on
loglog(matrix_size, LS_info_SVD(2, :), 'bd-', 'LineWidth', 2, 'MarkerFaceColor', 'blue'), hold on
loglog(matrix_size, LS_info_WSVD(3, :), 'ro:', 'LineWidth', 2), hold on
loglog(matrix_size, LS_info_SVD(3, :), 'bd:', 'LineWidth', 2), hold on
xlim([matrix_size(1)/1.2 1.2*matrix_size(end)])
xlabel('matrix size $3n$','FontSize',16,'Interpreter','latex'); 
ylabel('time(sec.)','FontSize',16);
legend('WSVD-CG', 'SVD-CG', 'construct WSVD', 'construct SVD', 'Location','southeast', 'FontSize', 14);
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
% saveas(gca, 'CG_time_FCC', 'epsc');

%% plot SC
load DATA/LS_info_SC.mat
n = size(LS_info_SVD, 2);
grid_nums = (1 : n)' * 10;
matrix_size = 3 * grid_nums.^3;
figure
loglog(matrix_size, LS_info_WSVD(2, :), 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'red'), hold on
loglog(matrix_size, LS_info_SVD(2, :), 'bd-', 'LineWidth', 2, 'MarkerFaceColor', 'blue'), hold on
loglog(matrix_size, LS_info_WSVD(3, :), 'ro:', 'LineWidth', 2), hold on
loglog(matrix_size, LS_info_SVD(3, :), 'bd:', 'LineWidth', 2), hold on
xlim([matrix_size(1)/1.2 1.2*matrix_size(end)])
xlabel('matrix size $3n$','FontSize',16,'Interpreter','latex'); 
ylabel('time(sec.)','FontSize',16);
legend('WSVD-CG', 'SVD-CG', 'construct WSVD', 'construct SVD', 'Location','southeast', 'FontSize', 14);
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
% saveas(gca, 'CG_time_SC', 'epsc');

%% plot HEX
load DATA/LS_info_HEX.mat
n = size(LS_info_SVD, 2);
grid_nums = (1 : n)' * 10;
matrix_size = 3 * grid_nums.^3;
figure
loglog(matrix_size, LS_info_WSVD(2, :), 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'red'), hold on
loglog(matrix_size, LS_info_SVD(2, :), 'bd-', 'LineWidth', 2, 'MarkerFaceColor', 'blue'), hold on
loglog(matrix_size, LS_info_WSVD(3, :), 'ro:', 'LineWidth', 2), hold on
loglog(matrix_size, LS_info_SVD(3, :), 'bd:', 'LineWidth', 2), hold on
xlim([matrix_size(1)/1.2 1.2*matrix_size(end)])
xlabel('matrix size $3n$','FontSize',16,'Interpreter','latex'); 
ylabel('time(sec.)','FontSize',16);
legend('WSVD-CG', 'SVD-CG', 'construct WSVD', 'construct SVD', 'Location','southeast', 'FontSize', 14);
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
% saveas(gca, 'CG_time_HEX', 'epsc');




