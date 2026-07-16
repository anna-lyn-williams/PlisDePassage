%% Plot distances
load('distancesaddlhfinal.mat', 'distancesaddlhfinal');

% fsaverage reference distances from the TOCS (mm), one per boundary.
% Uncomment the matching dotted-line overlay below to draw them.
tocs_to_ab0 = 23;
tocs_to_01  = 39;
tocs_to_12  = 56;
tocs_to_23  = 71;

f=figure('Name','Left Hemisphere Annectant Gyri Locations','color','w');

% Boundary 1
h1=histfit(distancesaddlhfinal(:,1), 5, 'kernel');
h1(2).Color = [0.8125 0.2422 0.4961];
delete(h1(1))
avg = mean(distancesaddlhfinal(:,1));
hold on;
line([avg avg], ylim, 'Color', [0.8125 0.2422 0.4961], 'LineWidth', 2, 'LineStyle', '--');   % Observation Group mean
%line([tocs_to_ab0 tocs_to_ab0], ylim, 'Color', [0.8125 0.2422 0.4961], 'LineWidth', 2, 'LineStyle', ':');   % fsaverage reference

% Boundary 2
h2=histfit(distancesaddlhfinal(:,2), 5, 'kernel');
h2(2).Color = [0.9453 0.6914 0.2382];
delete(h2(1))
avg = mean(distancesaddlhfinal(:,2));
hold on;
line([avg avg], ylim, 'Color', [0.9453 0.6914 0.2382], 'LineWidth', 2, 'LineStyle', '--');   % Observation Group mean
%line([tocs_to_01 tocs_to_01], ylim, 'Color', [0.9453 0.6914 0.2382], 'LineWidth', 2, 'LineStyle', ':');   % fsaverage reference

% Boundary 3
h3=histfit(distancesaddlhfinal(:,3), 5, 'kernel');
h3(2).Color = 'green';
delete(h3(1))
avg = mean(distancesaddlhfinal(:,3));
hold on;
line([avg avg], ylim, 'Color', 'green', 'LineWidth', 2, 'LineStyle', '--');   % Observation Group mean
%line([tocs_to_12 tocs_to_12], ylim, 'Color', 'green', 'LineWidth', 2, 'LineStyle', ':');   % fsaverage reference

% Boundary 4
h4=histfit(distancesaddlhfinal(:,4), 5, 'kernel');
h4(2).Color = [0.1133 0.5469 0.5781];
delete(h4(1))
avg = mean(distancesaddlhfinal(:,4));
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
std(distancesaddlhfinal)
mean(distancesaddlhfinal)

data = array2table(distancesaddlhfinal, ...
    'VariableNames', {'D1','D2','D3','D4'});
within = table( ...
    categorical({'D1'; 'D2'; 'D3'; 'D4'}), ...
    'VariableNames', {'DistanceCondition'});
rm = fitrm(data, 'D1-D4 ~ 1', ...
           'WithinDesign', within);
ranova(rm)
