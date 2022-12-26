%% plot iteration number of linear system
%% simple cubic lattice

clc 
clear 
close all

%% input

filename = 'DATA/LS_info_SC.txt'; 
T = readtable(filename);
data = table2array(T);

num = size(data, 1);
mesh = (20:10: 20 + 10*(num-1));

load DATA/ConditionNumbers_SC.mat

amg.iter = data(:,1); amg.ful_time = data(:,2); amg.CG_time = data(:,3);

Wsvd.iter = data(:,4); Wsvd.ful_time = data(:,5); Wsvd.pre = data(:,6); Wsvd.CG = Wsvd.ful_time - Wsvd.pre;

svd.iter = data(:,7); svd.ful_time = data(:,8); svd.pre = data(:,9); svd.CG = svd.ful_time - svd.pre;

amg.pre = amg.ful_time - amg.CG_time;

matrix_size = mesh;
% matrix_size = 3*mesh.^3;

%% Plot1: pre

% plot1: pre
semilogy(matrix_size, amg.pre,':s','Color','black',  'LineWidth',2);
hold on 
semilogy(matrix_size, Wsvd.pre,':o','Color','red', 'LineWidth',2);
hold on 
semilogy(matrix_size, svd.pre,':d','Color','blue', 'LineWidth',2);
hold on 

semilogy(matrix_size, amg.CG_time,'-s','Color','black', 'LineWidth',2, 'MarkerFaceColor', 'black');
hold on 
semilogy(matrix_size, Wsvd.CG,'-o','Color','red', 'LineWidth',2, 'MarkerFaceColor','red');
hold on 
semilogy(matrix_size, svd.CG,'-d','Color','blue', 'LineWidth',2, 'MarkerFaceColor','blue');
hold on 

xlabel('$n_1$','FontSize',16,'Interpreter','latex'); 
ylabel('time(sec.)','FontSize',14);
legend('construct AMG', 'construct WSVD', 'construct SVD', 'AMG-CG', 'WSVD-CG', 'SVD-CG', 'Location','southeast');

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

%% Plot2: iter

figure
% plot2: iter
clf reset
h1 = axes('Position', [0.12 0.12 0.82 0.82]);
plot(matrix_size, amg.iter,'--s','Color','black',  'LineWidth',2, 'MarkerFaceColor', 'black');
hold on 
plot(matrix_size, Wsvd.iter,'--o','Color','red', 'LineWidth',2, 'MarkerFaceColor','red');
hold on 
plot(matrix_size, svd.iter,'--d','Color','blue', 'LineWidth',2, 'MarkerFaceColor', 'blue');
hold on 

xlabel('$n_1$','FontSize',16,'Interpreter','latex'); 
ylabel('iteration number','FontSize',14);
legend('AMG-CG', 'WSVD-CG', 'SVD-CG', 'Location','northwest');

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

h2 = axes('Position', [0.65 0.3 0.25 0.25]);
plot(matrix_size, CD_SC, '+', 'Color', 'm', 'LineWidth', 2, 'MarkerFaceColor', 'm'), hold on
axis(h2, [10 170 min(CD_SC) - 10 max(CD_SC) + 10])





