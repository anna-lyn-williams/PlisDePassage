%% pRF eccentricity at annectant-gyrus boundaries
% Loads eccentricity at each boundary for the true data and the spin-test
% (control) null, summarises the mean per boundary, and compares peripheral
% (AB/0, 1/2) against foveal (0/1, 2/3) boundaries with a repeated-measures
% ANOVA.
%
% To run:
% 1) set trueDataDir and controlDataDir paths and 2) press 'Run'
% 
% Swap 'rh' for 'lh' in the filenames below to run the other
% hemisphere. Specifically, in load('rh.adults.eccentricity_table_AW.mat') and
% load('rh.adults.eccentricity_table_AW.mat') below.


clear; clc; close all;

%% Boundary columns, ordered AB/0, 0/1, 1/2, 2/3 (must match table variable names)
vars = {'A.Gyrus.ab0', 'A.Gyrus.01', 'A.Gyrus.12', 'A.Gyrus.23'};

%% True data

trueDataDir = 'path/to/annectant_gyri_stats';
cd(trueDataDir);
load('rh.adults.eccentricity_table_AW.mat');  % loads eccentricity_table
trueDataTable = eccentricity_table;

true_rows = 1:height(trueDataTable);

% Grab the relevant columns
true_data = trueDataTable{true_rows, vars};

% Coerce to numeric (tables may hold cells or doubles)
if iscell(true_data)
    true_numericData = NaN(size(true_data));
    for r = 1:size(true_data, 1)
        for c = 1:size(true_data, 2)
            if isempty(true_data{r, c})
                true_numericData(r, c) = NaN;
            else
                true_numericData(r, c) = true_data{r, c};
            end
        end
    end
else
    true_numericData = true_data;
end

% Keep only subjects with a value at all 4 boundaries
true_validRows   = all(~isnan(true_numericData), 2);
true_filteredData = true_numericData(true_validRows, :);

% Mean eccentricity at each boundary (AB/0, 0/1, 1/2, 2/3)
true_meanValues = mean(true_filteredData, 1);


%% Control data (spin-test null)

controlDataDir = 'path/to/spintest_stats';
cd(controlDataDir);
load('rh.adults.eccentricity_table_AW.mat');  % loads eccentricity_table
controlDataTable = eccentricity_table;

control_rows = 1:height(controlDataTable);

control_data = controlDataTable{control_rows, vars};

% Coerce to numeric (tables may hold cells or doubles)
if iscell(control_data)
    control_numericData = NaN(size(control_data));
    for r = 1:size(control_data, 1)
        for c = 1:size(control_data, 2)
            if isempty(control_data{r, c})
                control_numericData(r, c) = NaN;
            else
                control_numericData(r, c) = control_data{r, c};
            end
        end
    end
else
    control_numericData = control_data;
end

% Keep only subjects with a value at all 4 boundaries
control_validRows    = all(~isnan(control_numericData), 2);
control_filteredData = control_numericData(control_validRows, :);

% Mean eccentricity at each boundary
control_meanValues = mean(control_filteredData, 1);


% Optional: mean +/- SE across boundaries, true vs control (uncomment to draw)
true_SE    = std(true_filteredData, [], 1)    / sqrt(size(true_filteredData, 1));
control_SE = std(control_filteredData, [], 1) / sqrt(size(control_filteredData, 1));

xAxis = 1:size(true_filteredData, 2);  % 1:4 for AB/0, 0/1, 1/2, 2/3
xPatch = [1 2 3 4 4 3 2 1];

figure('Color', 'w'); hold on;

darkBlue = [17 28 255] / 255;   % HEX #111cff

% True data patch (mean +/- SE)
y_true = [true_meanValues + true_SE, fliplr(true_meanValues - true_SE)];
patch(xPatch, y_true, darkBlue, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(xAxis, true_meanValues, 'Color', darkBlue, 'LineWidth', 2, ...
     'DisplayName', 'True mean');

% Control data patch (mean +/- SE)
y_control = [control_meanValues + control_SE, fliplr(control_meanValues - control_SE)];
patch(xPatch, y_control, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(xAxis, control_meanValues, 'Color', 'k', 'LineWidth', 2, ...
     'DisplayName', 'Control mean');

xlabel('Annectant Gyrus Boundary', 'Color', 'k');
ylabel('Eccentricity', 'Color', 'k');
title('pRF Eccentricity at Anat Boundaries', 'Color', 'k');
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'AB/0', '0/1', '1/2', '2/3'},fontsize=20);

hold off;


%% Repeated-measures ANOVA: peripheral vs foveal boundaries
peripheral = mean(true_filteredData(:, [1 3]), 2); % AB0, 1/2
foveal     = mean(true_filteredData(:, [2 4]), 2); % 0/1, 2/3

nSubj = size(true_filteredData, 1);
subjects = (1:nSubj)';

% Build repeated-measures table
T = table(subjects, peripheral, foveal, ...
          'VariableNames', {'Subject', 'Peripheral', 'Foveal'});

% Within-subject factor
within = table({'Peripheral'; 'Foveal'}, 'VariableNames', {'Representation'});

% Fit repeated-measures model
rm = fitrm(T, 'Peripheral-Foveal ~ 1', 'WithinDesign', within);

% ANOVA for Representation
ranova_tbl = ranova(rm, 'WithinModel', 'Representation');

% Display results
disp(ranova_tbl)

%% Descriptive statistics
mean_peripheral = mean(peripheral)
sd_peripheral   = std(peripheral)

mean_foveal = mean(foveal)
sd_foveal   = std(foveal)
