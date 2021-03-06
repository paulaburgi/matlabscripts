% run cctest_avgvel_randBl.m for m2 and b2. 

% 
% h = figure; hold on; box on; 
% options = struct('x_axis', ccdates, 'error', 'std', 'handle', gcf, 'color_area', ...
%     [0.5 0.5 0.5], 'color_line', [0.5 0.5 0.5], 'alpha', 0.33, 'line_width', 1);
% yyaxis left
% plot_areaerrorbar(m2, options);
% yyaxis right
% hold on
% options = struct('x_axis', ccdates, 'error', 'std', 'handle', gcf, 'color_area', ...
%     [0 0.5 0], 'color_line', [0 0.5 0], 'alpha', 0.33, 'line_width', 1);
% plot_areaerrorbar(b2, options);
% ylim([-40 40])
% datetick; 
% 
% % for sentinel duration tests
% h = figure; hold on; box on; 
% options = struct('x_axis', d, 'error', 'std', 'handle', gcf, 'color_area', ...
%     [0.5 0.5 0.5], 'color_line', [0.5 0.5 0.5], 'alpha', 0.33, 'line_width', 1);
%     figure(options.handle);
% x_vector = [options.x_axis', fliplr(options.x_axis')];
% error1 = [ms_all-mm_all(:,1)]'; 
% error2 = [mm_all(:,2)-ms_all]';
% data_mean = ms_all';
% patch = fill(x_vector, [data_mean+error2, fliplr(data_mean-error1)], options.color_area);
% set(patch, 'edgecolor', [0.7 0.7 0.7]);
% set(patch, 'FaceAlpha', options.alpha);
% hold on;
% plot(options.x_axis, data_mean, '-', 'color', options.color_line, ...
%     'LineWidth', options.line_width);
% plot(options.x_axis, data_mean, '.', 'markersize', 7, 'color', options.color_line, ...
%     'LineWidth', options.line_width);
% d; 
% xlim([0 5.5]);


close all

% for sentinel linear vs non linear tests
h = figure; hold on; box on; 
yyaxis left
hold on; 
options = struct('x_axis', ccdates, 'error', 'std', 'handle', h, 'color_area', ...
    [0.5 0.5 0.5], 'color_line', [0.5 0.5 0.5], 'alpha', 0.2, 'line_width', 1, ...
    'markertype', 'o', 'markersize', 7);
    plot_areaerrorbar(m1, options); hold on; 
options = struct('x_axis', ccdates, 'error', 'std', 'handle', h, 'color_area', ...
    [0.8 0.8 0.8], 'color_line', [0.4 0.4 0.4], 'alpha', 0.2, 'line_width', 1, ...
    'markertype', 'x', 'markersize', 10);
    plot_areaerrorbar(m2, options); hold on; 
options = struct('x_axis', ccdates, 'error', 'std', 'handle', h, 'color_area', ...
    [0.3 0.3 0.3], 'color_line', [0.3 0.3 0.3], 'alpha', 0.2, 'line_width', 1, ...
    'markertype', 's', 'markersize', 5);
    plot_areaerrorbar(m3, options); hold on; 
ylim([-1 1])
yyaxis right
hold on
options = struct('x_axis', ccdates, 'error', 'std', 'handle', h, 'color_area', ...
    [0 0.5 0], 'color_line', [0 0.6 0], 'alpha', 0.2, 'line_width', 1, ...
    'markertype', 'o', 'markersize', 7);
    plot_areaerrorbar(b1, options); hold on; 
options = struct('x_axis', ccdates, 'error', 'std', 'handle', h, 'color_area', ...
    [0 0.8 0], 'color_line', [0 0.5 0], 'alpha', 0.2, 'line_width', 1, ...
    'markertype', 'x', 'markersize', 10);
    plot_areaerrorbar(b2, options); hold on; 
options = struct('x_axis', ccdates, 'error', 'std', 'handle', h, 'color_area', ...
    [0 0.3 0], 'color_line', [0 0.4 0], 'alpha', 0.2, 'line_width', 1, ...
    'markertype', 's', 'markersize', 5);
    plot_areaerrorbar(b3, options); hold on; 
datetick; 
ylim([-40 40])






