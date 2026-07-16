%% Plot distances
load('distancesaddrhfinal.mat', 'distancesaddrhfinal');

% fsaverage reference distances from the TOCS (mm), one per boundary.
% Uncomment the matching dotted-line overlay below to draw them.
tocs_to_ab0 = 21;
tocs_to_01  = 34;
tocs_to_12  = 45;
tocs_to_23  = 63;

f=figure('Name','Left Hemisphere Annectant Gyri Locations','color','w');

% Boundary 1
h1=histfit(distancesaddrhfinal(:,1), 5, 'kernel');
h1(2).Color = [0.8125 0.2422 0.4961];
delete(h1(1))
avg = mean(distancesaddrhfinal(:,1));
hold on;
line([avg avg], ylim, 'Color', [0.8125 0.2422 0.4961], 'LineWidth', 2, 'LineStyle', '--');   % Observation Group mean
%line([tocs_to_ab0 tocs_to_ab0], ylim, 'Color', [0.8125 0.2422 0.4961], 'LineWidth', 2, 'LineStyle', ':');   % fsaverage reference

% Boundary 2
h2=histfit(distancesaddrhfinal(:,2), 5, 'kernel');
h2(2).Color = [0.9453 0.6914 0.2382];
delete(h2(1))
avg = mean(distancesaddrhfinal(:,2));
hold on;
line([avg avg], ylim, 'Color', [0.9453 0.6914 0.2382], 'LineWidth', 2, 'LineStyle', '--');   % Observation Group mean
%line([tocs_to_01 tocs_to_01], ylim, 'Color', [0.9453 0.6914 0.2382], 'LineWidth', 2, 'LineStyle', ':');   % fsaverage reference

% Boundary 3
h3=histfit(distancesaddrhfinal(:,3), 5, 'kernel');
h3(2).Color = 'green';
delete(h3(1))
avg = mean(distancesaddrhfinal(:,3));
hold on;
line([avg avg], ylim, 'Color', 'green', 'LineWidth', 2, 'LineStyle', '--');   % Observation Group mean
%line([tocs_to_12 tocs_to_12], ylim, 'Color', 'green', 'LineWidth', 2, 'LineStyle', ':');   % fsaverage reference

% Boundary 4
h4=histfit(distancesaddrhfinal(:,4), 5, 'kernel');
h4(2).Color = [0.1133 0.5469 0.5781];
delete(h4(1))
avg = mean(distancesaddrhfinal(:,4));
hold on;
line([avg avg], ylim, 'Color', [0.1133 0.5469 0.5781], 'LineWidth', 2, 'LineStyle', '--');   % Observation Group mean
%line([tocs_to_23 tocs_to_23], ylim, 'Color', [0.1133 0.5469 0.5781], 'LineWidth', 2, 'LineStyle', ':');   % fsaverage reference

xticks('auto')
xlabel('Geodesic Distance from TOCS (mm)', 'FontSize', 28,'FontWeight', 'bold')
ylabel('Count', 'FontSize', 28,'FontWeight', 'bold')

set(gca, ...
    'XLim', [0 140], ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'LineWidth', 2, ...
    'FontSize', 22);


%% Standard deviations and repeated-measures ANOVA across boundaries
std(distancesaddrhfinal)
mean(distancesaddrhfinal)

data = array2table(distancesaddrhfinal, ...
    'VariableNames', {'D1','D2','D3','D4'});
within = table( ...
    categorical({'D1'; 'D2'; 'D3'; 'D4'}), ...
    'VariableNames', {'DistanceCondition'});
rm = fitrm(data, 'D1-D4 ~ 1', ...
           'WithinDesign', within);
ranova(rm)
