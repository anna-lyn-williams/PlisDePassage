function plotBoundaryDistances_rmanova(distanceType, subjectType, atlasComparison)
% Compare distances between the anatomical and probabilistic (Wang) atlas
% boundaries and test them with repeated-measures ANOVA.
%   distanceType    : 'centroid' or 'meanmin'
%   subjectType     : subject subset label (default 'all')
%   atlasComparison : atlas pair to compare (default 'anat_wang')
% Draws violin plots per hemisphere and a pooled LH+RH figure, and writes the
% ANOVA results to figures/stats.txt.
%
% Run with command: 
% plotBoundaryDistances_rmanova('centroid','all','anat_wang')

% Input handling
validTypes = {'centroid','meanmin'};
if ~ismember(lower(distanceType), validTypes)
    error('distanceType must be centroid or meanmin');
end
if nargin < 2 || isempty(subjectType),     subjectType    = 'all';       end
if nargin < 3 || isempty(atlasComparison), atlasComparison = 'anat_wang'; end

dataPath   = './';
outputPath = fullfile(dataPath,'figures');
if ~exist(outputPath,'dir'), mkdir(outputPath); end
addpath('./Violinplot-Matlab')

% Load data
datasetFiles = { ...
    'replication_adults_all_subjects_boundary_distances.csv',   'Replication (adults)'; ...
    'observation_adults_all_subjects_boundary_distances.csv',   'Observation'; ...
    'replication_children_all_subjects_boundary_distances.csv', 'Replication (children)'; ...
};

allTables = cell(size(datasetFiles,1), 1);
for d = 1:size(datasetFiles,1)
    tmp = readtable(fullfile(dataPath, datasetFiles{d,1}));
    tmp.dataset = repmat(datasetFiles(d,2), height(tmp), 1);
    if strcmpi(distanceType,'centroid')
        tmp.atlas1 = tmp.anat_centroid;
        tmp.atlas2 = tmp.wang_centroid;
    else
        tmp.atlas1 = tmp.anat_meanmin;
        tmp.atlas2 = tmp.wang_meanmin;
    end
    allTables{d} = tmp(:,{'subject','hemisphere','boundary','atlas1','atlas2','dataset'});
end

data = vertcat(allTables{:});
data.subjKey = cellstr(strcat(string(data.dataset),"__",string(data.subject)));

subjects    = unique(data.subjKey);
boundaries  = {'V3ab_IPS0','IPS0_IPS1','IPS1_IPS2','IPS2_IPS3'};
hemispheres = {'lh','rh'};

colors.atlas1 = [0.25 0.25 0.25];   % anat        -> dark gray
colors.atlas2 = [0.65 0.65 0.65];   % wang (prob) -> light gray

% Open stats file
statsFile = fullfile(outputPath,'stats.txt');
fileID    = fopen(statsFile,'w');

% Figure: lh & rh panels
fig_hemis = figure('Color','w','Position',[100 100 1000 500]);

atlas1_all = nan(numel(subjects), numel(boundaries), 2);
atlas2_all = nan(numel(subjects), numel(boundaries), 2);

for h = 1:2
    hemi = hemispheres{h};
    subplot(1,2,h); hold on

    [a1, a2] = extractData(data, subjects, boundaries, hemi);
    atlas1_all(:,:,h) = a1;
    atlas2_all(:,:,h) = a2;

    plotViolins(a1, a2, boundaries, colors);
    title(upper(hemi), 'FontSize', 24)

    runRMANOVA(a1, a2, subjects, boundaries, fileID, ...
        sprintf('%s Hemisphere', upper(hemi)));
end

exportgraphics(fig_hemis, fullfile(outputPath,'BoundaryDistances_hemis.pdf'), ...
    'ContentType','vector','BackgroundColor','none');

% Figure: both hemispheres combined
fig_avg = figure('Color','w','Position',[100 100 500 500]);
hold on

atlas1_pooled = [atlas1_all(:,:,1); atlas1_all(:,:,2)];
atlas2_pooled = [atlas2_all(:,:,1); atlas2_all(:,:,2)];

plotViolins(atlas1_pooled, atlas2_pooled, boundaries, colors);
title('LH + RH (Pooled)', 'FontSize', 24)

runRMANOVAWithHemi(atlas1_all, atlas2_all, subjects, boundaries, ...
    hemispheres, fileID, 'Combined LH+RH');

fclose(fileID);

exportgraphics(fig_avg, fullfile(outputPath,'BoundaryDistances_combined.pdf'), ...
    'ContentType','vector','BackgroundColor','none');

fprintf('Done. Stats saved to %s\n', statsFile);


% Local functions

function [atlas1_data, atlas2_data] = extractData(data, subjects, boundaries, hemi)

    atlas1_data = nan(numel(subjects), numel(boundaries));
    atlas2_data = nan(numel(subjects), numel(boundaries));

    for s = 1:numel(subjects)
        for b = 1:numel(boundaries)
            idx = strcmp(data.subjKey, subjects{s}) & ...
                  strcmp(data.hemisphere, hemi)      & ...
                  strcmp(data.boundary,   boundaries{b});
            if any(idx)
                atlas1_data(s,b) = data.atlas1(find(idx,1));
                atlas2_data(s,b) = data.atlas2(find(idx,1));
            end
        end
    end
end

% -------------------------------------------------------------------------
function plotViolins(a1, a2, boundaries, colors)
% Double half-violins matched to the compareGroupsSummary_bothhemis style:
% filled half-violin (kernel density) + small jittered points + mean line.
%   left half  = atlas1 (anat, dark gray)
%   right half = atlas2 (wang/prob, light gray)

    nB = size(a1,2);
    for b = 1:nB
        localHalfViolin(gca, b, a1(:,b), 'left',  colors.atlas1, 0.45);
        localHalfViolin(gca, b, a2(:,b), 'right', colors.atlas2, 0.45);
    end

    % mean trend lines across boundaries (kept from the original figure)
    x = 1:nB;
    plot(x, mean(a1,1,'omitnan'), 'LineWidth',3, 'Color',colors.atlas1)
    plot(x, mean(a2,1,'omitnan'), 'LineWidth',3, 'Color',colors.atlas2)

    set(gca,'FontSize',20,'LineWidth',3,'Box','off','TickLabelInterpreter','none')
    xlim([0.5, nB+0.5])
    xticks(1:nB)
    xticklabels(boundaries)
    ylabel('Distance (mm)')
end

% -------------------------------------------------------------------------
function runRMANOVA(a1, a2, subjects, boundaries, fileID, label)
% Single-hemisphere repeated measures ANOVA.
%
% Design:  within-subject factors = AtlasType (2) x Boundary (4)
%
% fitrm requires a wide table: one row per subject, one column per
% condition cell.  Column naming convention:
%   AT<atlasIdx>_B<boundaryIdx>   (e.g. AT1_B1 ... AT2_B4)
%
% Subjects with ANY missing cell are listwise-excluded so the design
% stays fully balanced (required by fitrm).

    fprintf(fileID, '\n========================================\n');
    fprintf(fileID, '%s - Repeated Measures ANOVA\n', label);
    fprintf(fileID, 'Within factors: AtlasType (2) x Boundary (4)\n');
    fprintf(fileID, '========================================\n');

    nS = numel(subjects);
    nB = numel(boundaries);

    % --- Build wide table (subjects x 8 condition columns) ---
    colNames = cell(1, 2*nB);
    dataWide = nan(nS, 2*nB);
    for b = 1:nB
        colNames{b}      = sprintf('AT1_B%d', b);   % atlas1 (anat)
        colNames{nB + b} = sprintf('AT2_B%d', b);   % atlas2 (wang)
        dataWide(:, b)      = a1(:, b);
        dataWide(:, nB + b) = a2(:, b);
    end

    wideT = array2table(dataWide, 'VariableNames', colNames);

    % Listwise deletion of subjects with any NaN
    complete = all(~isnan(dataWide), 2);
    wideT    = wideT(complete, :);
    nUsed    = sum(complete);

    fprintf(fileID, 'Subjects with complete data: %d / %d\n', nUsed, nS);

    if nUsed < 3
        fprintf(fileID, 'Insufficient subjects for rmanova — skipping.\n');
        return
    end

    % --- Within-design table ---
    % Two factors: AtlasType (1 or 2) and Boundary (1..4)
    factorAT = [ones(nB,1); 2*ones(nB,1)];
    factorB  = repmat((1:nB)', 2, 1);
    withinT  = table(categorical(factorAT), categorical(factorB), ...
        'VariableNames', {'AtlasType','Boundary'});

    % --- Fit and run ranova ---
    rm  = fitrm(wideT, [strjoin(colNames,',') ' ~ 1'], 'WithinDesign', withinT);
    rta = ranova(rm, 'WithinModel', 'AtlasType * Boundary');

    % --- Write results ---
    fprintf(fileID, '\n--- ranova Table ---\n');
    disp(rta)

    termNames = rta.Properties.RowNames;

    for i = 1:height(rta)
        fprintf(fileID, '%-40s  df=%.0f  F=%.3f  p=%.4g  pGG=%.4g\n', ...
            termNames{i}, rta.DF(i), rta.F(i), rta.pValue(i), rta.pValueGG(i));
    end

    % --- Post-hoc: AtlasType per Boundary if interaction is significant ---
    intIdx = contains(termNames, 'AtlasType') & contains(termNames, 'Boundary');
    if any(intIdx) && min(rta.pValue(intIdx)) < 0.05
        fprintf(fileID, '\n--- Post-hoc: AtlasType per Boundary (Bonferroni, interaction p<0.05) ---\n');
        runPostHocPerBoundary(a1(complete,:), a2(complete,:), ...
            subjects(complete), boundaries, fileID);
    end

    fprintf(fileID, '\n');
end

% -------------------------------------------------------------------------
function runRMANOVAWithHemi(atlas1_all, atlas2_all, subjects, boundaries, ...
                             hemispheres, fileID, label)
% Combined repeated measures ANOVA including Hemisphere as a third
% within-subject factor.
%
% Design:  AtlasType (2) x Boundary (4) x Hemisphere (2)  — fully crossed
% Column naming: AT<a>_B<b>_H<h>   (16 columns total)
%
% Subjects are excluded listwise if any of the 16 cells is NaN.

    fprintf(fileID, '\n========================================\n');
    fprintf(fileID, '%s - Repeated Measures ANOVA\n', label);
    fprintf(fileID, 'Within factors: AtlasType (2) x Boundary (4) x Hemisphere (2)\n');
    fprintf(fileID, '========================================\n');

    nS = numel(subjects);
    nB = numel(boundaries);
    nH = numel(hemispheres);
    nCols = 2 * nB * nH;   % 16

    colNames = cell(1, nCols);
    dataWide = nan(nS, nCols);

    col = 0;
    for a = 1:2
        for b = 1:nB
            for h = 1:nH
                col = col + 1;
                colNames{col} = sprintf('AT%d_B%d_H%d', a, b, h);
                if a == 1
                    dataWide(:, col) = atlas1_all(:, b, h);
                else
                    dataWide(:, col) = atlas2_all(:, b, h);
                end
            end
        end
    end

    wideT    = array2table(dataWide, 'VariableNames', colNames);
    complete = all(~isnan(dataWide), 2);
    wideT    = wideT(complete, :);
    nUsed    = sum(complete);

    fprintf(fileID, 'Subjects with complete data: %d / %d\n', nUsed, nS);

    if nUsed < 3
        fprintf(fileID, 'Insufficient subjects for rmanova — skipping.\n');
        return
    end

    % --- Within-design table (16 rows, one per column) ---
    factorAT   = [];  factorB = [];  factorH = [];
    for a = 1:2
        for b = 1:nB
            for h = 1:nH
                factorAT(end+1,1) = a;  %#ok<AGROW>
                factorB(end+1,1)  = b;  %#ok<AGROW>
                factorH(end+1,1)  = h;  %#ok<AGROW>
            end
        end
    end

    withinT = table(categorical(factorAT), categorical(factorB), categorical(factorH), ...
        'VariableNames', {'AtlasType','Boundary','Hemisphere'});

    % --- Fit and run ranova ---
    rm  = fitrm(wideT, [strjoin(colNames,',') ' ~ 1'], 'WithinDesign', withinT);
    rta = ranova(rm, 'WithinModel', 'AtlasType * Boundary * Hemisphere');

    % --- Write results ---
    fprintf(fileID, '\n--- ranova Table ---\n');
    disp(rta)

    termNames = rta.Properties.RowNames;
    for i = 1:height(rta)
        fprintf(fileID, '%-50s  df=%.0f  F=%.3f  p=%.4g  pGG=%.4g\n', ...
            termNames{i}, rta.DF(i), rta.F(i), rta.pValue(i), rta.pValueGG(i));
    end

    % --- Post-hoc: AtlasType per Boundary (pooling hemispheres) if interaction significant ---
    intIdx = contains(termNames, 'AtlasType') & contains(termNames, 'Boundary') ...
           & ~contains(termNames, 'Hemisphere');
    if any(intIdx) && min(rta.pValue(intIdx)) < 0.05
        fprintf(fileID, '\n--- Post-hoc: AtlasType per Boundary (Bonferroni, interaction p<0.05) ---\n');
        a1_pool = [atlas1_all(:,:,1); atlas1_all(:,:,2)];
        a2_pool = [atlas2_all(:,:,1); atlas2_all(:,:,2)];
        subj_pool = [subjects; subjects];
        complete_pool = all(~isnan([a1_pool, a2_pool]), 2);
        runPostHocPerBoundary(a1_pool(complete_pool,:), a2_pool(complete_pool,:), ...
            subj_pool(complete_pool), boundaries, fileID);
    end

    fprintf(fileID, '\n');
end

% -------------------------------------------------------------------------
function runPostHocPerBoundary(a1, a2, subjects, boundaries, fileID)
% For each boundary: paired t-test across subjects (anat vs wang).
% Bonferroni correction across the 4 boundaries.
%
% fitrm/ranova would be overkill here — a paired t-test per boundary
% is the standard simple follow-up when you only have 2 atlas levels.

    nB    = numel(boundaries);
    alpha = 0.05 / nB;
    fprintf(fileID, '  (Bonferroni-corrected alpha = %.4f for %d comparisons)\n', alpha, nB);

    for b = 1:nB
        d = a1(:,b) - a2(:,b);           % anat minus wang, per subject
        d = d(~isnan(d));                 % exclude subjects missing this boundary

        if numel(d) < 3
            fprintf(fileID, '  %-18s  too few observations\n', boundaries{b});
            continue
        end

        [~, p, ~, stats] = ttest(d);
        m1  = mean(a1(:,b),'omitnan');
        m2  = mean(a2(:,b),'omitnan');
        sig = '';
        if p < alpha,  sig = ' *';  end
        if p < 0.001,  sig = ' **'; end

        fprintf(fileID, '  %-18s  anat=%.3f  wang=%.3f  t(%d)=%.3f  p=%.4g%s\n', ...
            boundaries{b}, m1, m2, stats.df, stats.tstat, p, sig);
    end
end

% -------------------------------------------------------------------------
function localHalfViolin(ax, xCenter, vals, side, faceColor, faceAlpha)
% One half violin (kernel density) + small jittered points + mean line.
% Matches the style used in compareGroupsSummary_bothhemis.m.
    maxWidth = 0.35;                 % max horizontal extent of the half violin
    if strcmpi(side,'left'), s = -1; else, s = 1; end

    vals = vals(:);
    vals = vals(~isnan(vals));
    if isempty(vals), return; end

    % --- density patch (needs >=2 distinct values) ---
    if numel(vals) >= 2 && (max(vals) > min(vals))
        [f, xi] = ksdensity(vals);
        f  = f / max(f) * maxWidth;
        px = xCenter + s * f(:);
        polyX = [px;    repmat(xCenter, numel(xi), 1)];
        polyY = [xi(:); flipud(xi(:))];
        patch(ax, polyX, polyY, faceColor, ...
            'FaceAlpha', faceAlpha, 'EdgeColor', faceColor, 'LineWidth', 1.25);
    end

    % --- individual subject points, jittered onto the same side ---
    jit = s * (0.03 + 0.06 * rand(numel(vals),1));
    scatter(ax, xCenter + jit, vals, 14, faceColor, 'filled', ...
        'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor','none');

    % --- mean line ---
    mval = mean(vals);
    if s < 0
        plot(ax, [xCenter - maxWidth, xCenter], [mval mval], ...
            'Color', faceColor, 'LineWidth', 3);
    else
        plot(ax, [xCenter, xCenter + maxWidth], [mval mval], ...
            'Color', faceColor, 'LineWidth', 3);
    end
end

end  % end main function