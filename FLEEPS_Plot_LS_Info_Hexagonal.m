%% plot iteration number of linear system
%% hexagonal lattice

clc 
clear 
close all

%% input

filename = 'DATA/LS_info_Hexagonal.txt'; 
T = readtable(filename);
data = table2array(T);

num = size(data, 1);
mesh = (20:10: 20 + 10*(num-1));

load DATA/ConditionNumbers_Hexagonal.mat

agmg.iter = data(:,1); agmg.CG_time= data(:,2); agmg.pre = data(:,3);

Wsvd.iter = data(:,4); Wsvd.CG = data(:,5); Wsvd.pre = data(:,6); 

svd.iter = data(:,7); svd.CG = data(:,8); svd.pre = data(:,9); 

matrix_size = mesh;
% matrix_size = 3*mesh.^3;

%% Plot1: pre

% plot1: pre
semilogy(matrix_size, agmg.pre,':s','Color','black',  'LineWidth',2);
hold on 
semilogy(matrix_size, Wsvd.pre,':o','Color','red', 'LineWidth',2);
hold on 
semilogy(matrix_size, svd.pre,':d','Color','blue', 'LineWidth',2);
hold on 

semilogy(matrix_size, agmg.CG_time,'-s','Color','black', 'LineWidth',2, 'MarkerFaceColor', 'black');
hold on 
semilogy(matrix_size, Wsvd.CG,'-o','Color','red', 'LineWidth',2, 'MarkerFaceColor','red');
hold on 
semilogy(matrix_size, svd.CG,'-d','Color','blue', 'LineWidth',2, 'MarkerFaceColor','blue');
hold on 

xlabel('$n_1$','FontSize',16,'Interpreter','latex'); 
ylabel('time(sec.)','FontSize',14);
legend('construct AGMG', 'construct WSVD', 'construct SVD', 'AGMG-CG', 'WSVD-CG', 'SVD-CG', 'Location','southeast');

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
plot(matrix_size, agmg.iter,'--s','Color','black',  'LineWidth',2, 'MarkerFaceColor', 'black');
hold on 
plot(matrix_size, Wsvd.iter,'--o','Color','red', 'LineWidth',2, 'MarkerFaceColor','red');
hold on 
plot(matrix_size, svd.iter,'--d','Color','blue', 'LineWidth',2, 'MarkerFaceColor', 'blue');
hold on 

xlabel('$n_1$','FontSize',16,'Interpreter','latex'); 
ylabel('iteration number','FontSize',14);
legend('AGMG-CG', 'WSVD-CG', 'SVD-CG', 'Location','northwest');

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
plot(matrix_size, CD_HEX, '+', 'Color', 'm', 'LineWidth', 2, 'MarkerFaceColor', 'm'), hold on
axis(h2, [10 170 min(CD_HEX) - 10 max(CD_HEX) + 10])





