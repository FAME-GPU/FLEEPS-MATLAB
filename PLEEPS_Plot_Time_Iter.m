clear
clc
close all

lattice = 'FCC'; %'FCC', 'SC', 'HEX'
switch lattice
    case 'FCC'
        path_string = 'GXWKGLUWLK';
    case 'HEX'
        path_string = 'GMKGALHA';
    case 'SC'
        path_string = 'GXMGRX';
end
part_num = 10;

load(['DATA/Comput_Info_' lattice '_SVD.mat'])
cpu_time_SVD = comput_info.cpu_time';
wave_vec_num = length(cpu_time_SVD);
es_iter_SVD = zeros(wave_vec_num, 1);
ls_iter_SVD = zeros(wave_vec_num, 1);
for i = 1 : wave_vec_num
    tmp = comput_info.LS_iter{1, i};
    tmp = tmp(find(tmp > 0)); 
    es_iter_SVD(i) = length(tmp);
    if es_iter_SVD(i) > 200
        % Compute the actual eigs iteration for the expicit deflation for Gamma point 
        es_iter_SVD(i) = es_iter_SVD(i) / 4;
    end
    ls_iter_SVD(i) = mean(tmp);
end

load(['DATA/Comput_Info_' lattice '_WSVD.mat'])
cpu_time_WSVD = comput_info.cpu_time';
wave_vec_num = length(cpu_time_WSVD);
es_iter_WSVD = zeros(wave_vec_num, 1);
ls_iter_WSVD = zeros(wave_vec_num, 1);
for i = 1 : wave_vec_num
    tmp = comput_info.LS_iter{1, i};
    tmp = tmp(find(tmp > 0)); 
    es_iter_WSVD(i) = length(tmp);
    if es_iter_WSVD(i) > 200
        % Compute the actual eigs iteration for the expicit deflation for Gamma point 
        es_iter_WSVD(i) = es_iter_WSVD(i) / 4;
    end
    ls_iter_WSVD(i) = mean(tmp);
end

%% plot
hax = axes(figure);
part_len = 0.2; % Set width of partition bar
% Set XData
part_bar_idx = [];
vertex_idx   = 1;
XData        = 0;
part_bar_num = 0;
label_symbol_array{1} = path_string(1);
if path_string(1) == '|'
    error('Invalid path string format');
end
for i = 2:length(path_string)
    if (path_string(i-1) ~= '|') && (path_string(i) ~= '|')
        vertex_idx = [ vertex_idx, (vertex_idx(end)+part_num-1)];
        tmp = linspace( XData(end), (XData(end)+1), part_num);
        XData = [XData, tmp(2:end)];
        label_symbol_array{end+1} = path_string(i);
    elseif (path_string(i-1) == '|') && (path_string(i) ~= '|')
        vertex_idx = [ vertex_idx, (vertex_idx(end)+1) ];
        label_symbol_array{end+1} = path_string(i);
    elseif (path_string(i) == '|')
        part_bar_num = part_bar_num + 1;
        part_bar_idx = [part_bar_idx, [vertex_idx(end);vertex_idx(end)+1]];
        XData = [XData, (XData(end)+part_len)];
    end
end
for i = 1 : length(label_symbol_array)
    if label_symbol_array{i} == 'G'
        label_symbol_array{i} = 'Γ';
    end
end

% figure
% eigs_time
hax1 = subplot(3, 1, 1);
hold(hax1, 'on');
YData1 = cpu_time_WSVD / 1000;
YData2 = cpu_time_SVD / 1000;
plot(hax1, XData, YData1(1 : length(XData)), 'ro', 'MarkerSize', 5, 'LineWidth', 1.0), hold on
plot(hax1, XData, YData2(1 : length(XData)), 'bd', 'MarkerSize', 5, 'LineWidth', 1.0);
set(hax1, 'FontSize', 14);
set(hax1, 'XTick', XData(vertex_idx));
set(hax1, 'XTickLabel', upper(label_symbol_array));  
axis(hax1, [min(XData), max(XData), 0.6*max(mink(YData1, 3)), 1.2*min(maxk(YData2, 3))]);
grid(hax1, 'on');
% ytickformat('%d')
ylabel('time(10^3 sec)');

% eigs_iter
hax2 = subplot(3, 1, 2);
hold(hax2, 'on');
YData1 = es_iter_WSVD;
YData2 = es_iter_SVD;
plot(hax2, XData, YData1(1 : length(XData)), 'ro', 'MarkerSize', 5, 'LineWidth', 1.0), hold on
plot(hax2, XData, YData2(1 : length(XData)), 'bd', 'MarkerSize', 5, 'LineWidth', 1.0);
set(hax2, 'FontSize', 14);
set(hax2, 'XTick', XData(vertex_idx));
set(hax2, 'XTickLabel', upper(label_symbol_array));  
axis(hax2, [min(XData), max(XData), 0.6*min(YData1), 1.2*max(YData2)]);
grid(hax2, 'on');
ylabel('#iter(eigs)');

% ls_iter
hax3 = subplot(3, 1, 3);
hold(hax3, 'on');
YData1 = ls_iter_WSVD;
YData2 = ls_iter_SVD;
plot(hax3, XData, YData1(1 : length(XData)), 'ro', 'MarkerSize', 5, 'LineWidth', 1.0), hold on
plot(hax3, XData, YData2(1 : length(XData)), 'bd', 'MarkerSize', 5, 'LineWidth', 1.0);
set(hax3, 'FontSize', 14);
set(hax3, 'XTick', XData(vertex_idx));
set(hax3, 'XTickLabel', upper(label_symbol_array));  
axis(hax3, [min(XData), max(XData), 0.6*max(mink(YData1, 3)), 1.2*min(maxk(YData2, 3))]);
grid(hax3, 'on');
% ytickformat('%d')
ylabel('#iter(pcg)');
xlabel('wave vector \bf{k}');
legend('WSVD', 'SVD', 'Location', 'best', 'NumColumns', 2);
