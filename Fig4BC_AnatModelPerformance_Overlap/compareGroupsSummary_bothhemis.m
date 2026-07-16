function resultsTable = compareGroupsSummary_bothhemis(group1, group2, cohort, ageGroup, varargin)
% Compare two ROI-labelling methods (group1 vs group2) across dorsal-stream
% ROIs. Hemispheres are pooled by averaging LH and RH per subject per ROI:
%   X = (X_lh + X_rh) / 2
% One figure is produced per metric set (Hits/Misses/FA and Dice/d'/Jaccard),
% pairing a per-ROI line plot with distribution violins.
%
% Expected data files: <cohort>.<hemisphere>.<ageGroup>.<group>.<Metric>.mat
%
% The value stored as Accuracy.mat is the Jaccard index,
%   Jaccard = TP / (TP + FP + FN); it is loaded, shown and reported as Jaccard.
%
% Returns resultsTable (one row per metric, Hemisphere = 'both'). For each
% within-subject effect (Group, ROI, Group x ROI) it reports the p-value, F
% statistic, and numerator/denominator degrees of freedom, followed by the
% per-ROI t-test p-values. F(df1, df2) and p are also printed to the console.
%
% Commands to run: 
% compareGroupsSummary_bothhemis('Anat','Wang','Observation','adults')
% compareGroupsSummary_bothhemis('Anat','Wang','ReplicationCombined','combined')


%% === Input Validation ===
validGroups    = {'Anat', 'Wang'};
validCohorts   = {'Observation', 'Replication', 'ReplicationCombined'};
validAgeGroups = {'adults', 'children', 'combined'};

if ~ismember(group1, validGroups) || ~ismember(group2, validGroups)
    error('Invalid group. Choose from: Anat, anat_funcPrior, Wang');
end
if ~ismember(cohort, validCohorts)
    error('Invalid cohort. Choose from: Observation, Replication, ReplicationCombined');
end
if ~ismember(ageGroup, validAgeGroups)
    error('Invalid age group. Choose from: adults, children, combined');
end

% enforce combined pairing rules
if strcmp(cohort,'ReplicationCombined') && ~strcmp(ageGroup,'combined')
    error('For cohort ReplicationCombined, set ageGroup = combined');
end
if ~strcmp(cohort,'ReplicationCombined') && strcmp(ageGroup,'combined')
    error('ageGroup = combined is only valid for cohort ReplicationCombined');
end

%% === Parse Name–Value Parameters ===
defaultBaseDir = '.'; %_funcConstrained
p = inputParser;
addParameter(p, 'baseDir',  defaultBaseDir, @(x) ischar(x) || isstring(x));
addParameter(p, 'runStats', true,           @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});

baseDir  = char(p.Results.baseDir);
runStats = logical(p.Results.runStats);

%% === Setup ===
% ONLY PLOT IPS0–IPS3
variablesToPlot = {'IPS0', 'IPS1', 'IPS2', 'IPS3'};

% subject counts (define once)
n_Replication_adults   = 21;
n_Observation_adults   = 33;
n_Replication_children = 21;

% Filename pattern: <cohort>.<hemisphere>.<ageGroup>.<group>.<Metric>.mat
fileNames = struct();
fileNames.Dice                 = struct('fileFmt', '%s.%s.%s.%s.diceCoeficients.mat',     'VarName', 'diceTable');
fileNames.DPrime               = struct('fileFmt', '%s.%s.%s.%s.DPrime.mat',               'VarName', 'dPrimeTable');

% Stored under the "Accuracy" key but holds the Jaccard index
fileNames.Accuracy             = struct('fileFmt', '%s.%s.%s.%s.Accuracy.mat',             'VarName', 'accuracyTable');

fileNames.Hits                 = struct('fileFmt', '%s.%s.%s.%s.Hits.mat',                 'VarName', 'hitsTable');
fileNames.Misses               = struct('fileFmt', '%s.%s.%s.%s.Misses.mat',               'VarName', 'missesTable');

% --- split false alarms ---
fileNames.FalseAlarmsConfusion = struct('fileFmt', '%s.%s.%s.%s.FalseAlarmsConfusion.mat', 'VarName', 'falseAlarmsConfusionTable');
fileNames.FalseAlarmsBlank     = struct('fileFmt', '%s.%s.%s.%s.FalseAlarmsBlank.mat',     'VarName', 'falseAlarmsBlankTable');
fileNames.FalseAlarms          = struct('fileFmt', '%s.%s.%s.%s.FalseAlarms.mat',          'VarName', 'falseAlarmsTable');
% (FalseAlarmsTotal is computed, no file)

% paths per metric
paths = struct();
paths.Dice                 = fullfile(baseDir, 'Dice');
paths.DPrime               = fullfile(baseDir, 'DPrime');
paths.Accuracy             = fullfile(baseDir, 'Accuracy'); % Jaccard stored here
paths.Hits                 = fullfile(baseDir, 'HitsMissesFalseAlarms');
paths.Misses               = fullfile(baseDir, 'HitsMissesFalseAlarms');
paths.FalseAlarmsConfusion = fullfile(baseDir, 'HitsMissesFalseAlarms');
paths.FalseAlarmsBlank     = fullfile(baseDir, 'HitsMissesFalseAlarms');
paths.FalseAlarms          = fullfile(baseDir, 'HitsMissesFalseAlarms');

% Two metric sets, one figure each: Hits/Misses/FalseAlarms and Dice/d'/Jaccard
metricGroups = {{'Hits','Misses','FalseAlarms'}, {'Dice','DPrime','Accuracy'}};
metricLabels = {{'Hits','Misses','False Alarms'}, {'Dice','d''','Jaccard'}};

% line styles: group1 solid, group2 dotted
lineStyles = {'-', ':'};

% --- Metric colors (HEX requested) ---
% jaccard and false alarms: #49006a
% dice and hits:           #ae017e
% d' and misses:           #f768a1
metricHex = containers.Map();
metricHex('Accuracy')    = '#49006a'; % Jaccard (stored as Accuracy)
metricHex('FalseAlarms') = '#49006a';

metricHex('Dice')        = '#ae017e';
metricHex('Hits')        = '#ae017e';

metricHex('DPrime')      = '#f768a1';
metricHex('Misses')      = '#f768a1';

% results: HemLabel | Metric | p_Group | p_ROI | p_Group_x_ROI | p_ttests(vec)
resultsSummary = {};

%% === Helper: legend-friendly group names ===
    function nm = groupPrettyName(g)
        g = char(g);
        if strcmpi(g,'Wang')
            nm = 'prob';
        elseif strcmpi(g,'Anat')
            nm = 'anat';
        elseif strcmpi(g,'anat_funcPrior')
            nm = 'anat_funcPrior';
        else
            nm = g;
        end
    end

prettyG1 = groupPrettyName(group1);
prettyG2 = groupPrettyName(group2);

%% === Helper: HEX -> RGB (0..1) ===
    function rgb = hex2rgb(hex)
        hex = char(hex);
        if startsWith(hex,'#'), hex = hex(2:end); end
        if numel(hex) ~= 6
            error('hex2rgb: hex must be 6 characters (e.g., #49006a).');
        end
        rgb = [hex2dec(hex(1:2)), hex2dec(hex(3:4)), hex2dec(hex(5:6))] / 255;
    end

%% === Helper: robust table->numeric conversion (subjects x ROI), selecting ROI cols by name ===
    function Y = tableToNumeric(tab, n, roiNames, strict)
        % strict (default true): error if a requested ROI column is absent.
        % strict = false: missing ROI columns are filled with NaN. This is
        % used by the violin panel, which averages across whatever regions
        % are present (e.g. V3ab may not exist in every table).
        if nargin < 4 || isempty(strict), strict = true; end

        roiNames = cellstr(roiNames);
        vnames   = tab.Properties.VariableNames;

        Y = nan(n, numel(roiNames));

        for k = 1:numel(roiNames)
            ii = find(strcmp(vnames, roiNames{k}), 1);
            if isempty(ii)
                if strict
                    error('Missing ROI column "%s" in table. Found columns: %s', ...
                        roiNames{k}, strjoin(vnames, ', '));
                else
                    warning('tableToNumeric:missingROI', ...
                        'ROI column "%s" not found; filling with NaN.', roiNames{k});
                    continue;   % leave this column as NaN
                end
            end

            C = table2cell(tab(1:n, ii));   % n x 1 cell for this ROI
            for r = 1:n
                v = C{r};
                if isempty(v)
                    Y(r,k) = NaN;
                elseif isnumeric(v) && isscalar(v)
                    Y(r,k) = v;
                elseif iscell(v) && numel(v)==1 && isnumeric(v{1}) && isscalar(v{1})
                    Y(r,k) = v{1};
                elseif ischar(v) || isstring(v)
                    Y(r,k) = str2double(v);   % NaN if non-numeric
                else
                    Y(r,k) = NaN;
                end
            end
        end
    end

%% === Helper: load numeric matrix for a single cohort/ageGroup ===
    function Y = loadOne(metric, hemisphere, group, coh, ag, n, roiNames, strictCols)
        if nargin < 7 || isempty(roiNames),   roiNames   = variablesToPlot; end
        if nargin < 8 || isempty(strictCols), strictCols = true;            end

        if strcmp(metric, 'FalseAlarmsTotal')
            Yc = loadOne('FalseAlarmsConfusion', hemisphere, group, coh, ag, n, roiNames, strictCols);
            Yb = loadOne('FalseAlarmsBlank',     hemisphere, group, coh, ag, n, roiNames, strictCols);
            Y  = Yc + Yb;
            return;
        end

        dataPath = paths.(metric);
        fileFmt  = fileNames.(metric).fileFmt;
        varName  = fileNames.(metric).VarName;

        f = fullfile(dataPath, sprintf(fileFmt, coh, hemisphere, ag, group));
        if ~exist(f, 'file')
            error('Missing file:\n%s\n\n(Expecting pattern: <cohort>.<hemi>.<ageGroup>.<group>.<metric>.mat)', f);
        end

        tmp = load(f);
        if ~isfield(tmp, varName)
            error('File exists but missing variable "%s":\n%s', varName, f);
        end
        tab = tmp.(varName);

        Y = tableToNumeric(tab, n, roiNames, strictCols);
    end

%% === Helper: load data for the requested pool (normal or combined) ===
    function X = loadMetricForPool(metric, hemisphere, group, roiNames, strictCols)
        if nargin < 4 || isempty(roiNames),   roiNames   = variablesToPlot; end
        if nargin < 5 || isempty(strictCols), strictCols = true;            end

        if strcmp(cohort, 'ReplicationCombined')
            X = [ ...
                loadOne(metric, hemisphere, group, 'Replication', 'adults',   n_Replication_adults,   roiNames, strictCols); ...
                loadOne(metric, hemisphere, group, 'Replication', 'children', n_Replication_children, roiNames, strictCols) ...
            ];
        else
            if strcmp(cohort,'Observation') && strcmp(ageGroup,'adults')
                n = n_Observation_adults;
            elseif strcmp(cohort,'Replication') && strcmp(ageGroup,'adults')
                n = n_Replication_adults;
            elseif strcmp(cohort,'Replication') && strcmp(ageGroup,'children')
                n = n_Replication_children;
            else
                error('Invalid cohort/ageGroup combination for subject count');
            end
            X = loadOne(metric, hemisphere, group, cohort, ageGroup, n, roiNames, strictCols);
        end
    end

%% === Helper: hemisphere-averaged data (LH/RH averaged per subject per ROI) ===
    function X = loadMetricHemiAveraged(metric, group, roiNames, strictCols)
        if nargin < 3 || isempty(roiNames),   roiNames   = variablesToPlot; end
        if nargin < 4 || isempty(strictCols), strictCols = true;            end
        X_lh = loadMetricForPool(metric, 'lh', group, roiNames, strictCols);
        X_rh = loadMetricForPool(metric, 'rh', group, roiNames, strictCols);

        if size(X_lh,1) ~= size(X_rh,1)
            error('LH and RH subject counts differ for %s (%s).', metric, group);
        end
        if size(X_lh,2) ~= size(X_rh,2)
            error('LH and RH ROI counts differ for %s (%s).', metric, group);
        end

        X = (X_lh + X_rh) ./ 2;
    end

%% === VIOLIN PANEL (subplot 2) CONFIG ===
% Regions averaged (per subject) for each half-violin. V3ab is included when
% present in the data table; tables that lack it simply skip that column
% (strictCols = false) and the average falls back to the regions that exist.
violinRegions = [{'V3ab'}, variablesToPlot];   % V3ab, IPS0, IPS1, IPS2, IPS3

% Left->right order of the three violins in subplot 2, plus their tick labels.
%   metricGroups{1} = Hits / Misses / FalseAlarms
%   metricGroups{2} = Dice / DPrime / Accuracy(=Jaccard)
% For the Dice/d'/Jaccard figure the requested order is Jaccard, Dice, d'.
violinOrder      = { {'Hits','Misses','FalseAlarms'}, {'Accuracy','Dice','DPrime'} };
violinTickLabels = { {'Hits','Misses','FA'},          {'Jaccard','Dice','d'''}     };
violinXPositions = [1 2 3];

%% === MAIN LOOP: ONE figure per metric set (LH/RH averaged) ===
%  subplot 1 = ROI line plot (original) ; subplot 2 = double violins
for plotIdx = 1:2
    figure('Color','w','Position',[100 100 1200 600]);

    % Manual flush axes so the two panels read as ONE continuous plot:
    %   left  panel (line plot) = 2/3 of the combined width
    %   right panel (violins)   = 1/3 of the combined width
    %   right panel starts exactly where the left ends -> no whitespace between.
    axPosL = [0.11 0.16 0.56 0.74];   % [left bottom width height]
    axPosR = [0.67 0.16 0.28 0.74];   % left = 0.11+0.56 (flush); width = 0.56/2

    % ================= SUBPLOT 1: ROI line plot (original) =================
    ax1 = axes('Position', axPosL);
    hold(ax1,'on');

    legendHandles = [];
    legendLabels  = {};

    for m = 1:numel(metricGroups{plotIdx})
        metric = metricGroups{plotIdx}{m};
        label  = metricLabels{plotIdx}{m};

        % display "Accuracy" as "Jaccard" since Jaccard Index is actually
        % what is computed.
        if strcmp(metric, 'Accuracy')
            label = 'Jaccard';
        end

        g1 = loadMetricHemiAveraged(metric, group1);
        g2 = loadMetricHemiAveraged(metric, group2);

        if ~isKey(metricHex, metric)
            error('No color defined for metric "%s". Add it to metricHex.', metric);
        end
        thisColor = hex2rgb(metricHex(metric));

        h1 = plot(mean(g1,1,'omitnan'), 'LineStyle', lineStyles{1}, 'Color', thisColor, 'LineWidth', 2);
        h2 = plot(mean(g2,1,'omitnan'), 'LineStyle', lineStyles{2}, 'Color', thisColor, 'LineWidth', 2);

        legendHandles = [legendHandles h1 h2]; 
        legendLabels  = [legendLabels {sprintf('%s - %s', label, prettyG1), ...
                                       sprintf('%s - %s', label, prettyG2)}]; 

        if runStats
            try
                n_subj = size(g1, 1); 
                
                % Combined data: each subject has 8 measurements (4 ROIs x 2 methods)
                T = array2table([g1, g2], 'VariableNames', ...
                    [strcat(variablesToPlot, '_g1'), strcat(variablesToPlot, '_g2')]);
                
                % Within design: ROI x Group (both within)
                withinROI    = [repmat(variablesToPlot', 2, 1)];
                withinGroup  = [repmat({'g1'}, 4, 1); repmat({'g2'}, 4, 1)];
                within = table(categorical(withinROI), categorical(withinGroup), ...
                    'VariableNames', {'ROI', 'Group'});
                
                rm = fitrm(T, 'IPS0_g1-IPS3_g2 ~ 1', 'WithinDesign', within);
                ranovatbl = ranova(rm, 'WithinModel', 'ROI*Group');
                rows = string(ranovatbl.Row);   % needed by the idx_* lookups below
                                
                

                p_group       = NaN;  F_group = NaN;  df1_group = NaN;  df2_group = NaN;
                p_roi         = NaN;  F_roi   = NaN;  df1_roi   = NaN;  df2_roi   = NaN;
                p_interaction = NaN;  F_inter = NaN;  df1_inter = NaN;  df2_inter = NaN;

                % For each effect: F and p come from the effect row; the
                % numerator df is that row's DF and the denominator df is the
                % DF of its matching Error(...) row.
                idx_group     = find(contains(rows,"Group") & ~contains(rows,"ROI") & ~contains(rows,"Error"), 1);
                idx_group_err = find(contains(rows,"Group") & ~contains(rows,"ROI") &  contains(rows,"Error"), 1);
                if ~isempty(idx_group)
                    p_group   = ranovatbl.pValue(idx_group);
                    F_group   = ranovatbl.F(idx_group);
                    df1_group = ranovatbl.DF(idx_group);
                    if ~isempty(idx_group_err), df2_group = ranovatbl.DF(idx_group_err); end
                end

                idx_roi     = find(contains(rows,"ROI") & ~contains(rows,"Group") & ~contains(rows,"Error") & ~contains(rows,"(Intercept)"), 1);
                idx_roi_err = find(contains(rows,"ROI") & ~contains(rows,"Group") &  contains(rows,"Error"), 1);
                if ~isempty(idx_roi)
                    p_roi   = ranovatbl.pValue(idx_roi);
                    F_roi   = ranovatbl.F(idx_roi);
                    df1_roi = ranovatbl.DF(idx_roi);
                    if ~isempty(idx_roi_err), df2_roi = ranovatbl.DF(idx_roi_err); end
                end

                idx_inter     = find(contains(rows,"ROI") & contains(rows,"Group") & ~contains(rows,"Error"), 1);
                idx_inter_err = find(contains(rows,"ROI") & contains(rows,"Group") &  contains(rows,"Error"), 1);
                if ~isempty(idx_inter)
                    p_interaction = ranovatbl.pValue(idx_inter);
                    F_inter       = ranovatbl.F(idx_inter);
                    df1_inter     = ranovatbl.DF(idx_inter);
                    if ~isempty(idx_inter_err), df2_inter = ranovatbl.DF(idx_inter_err); end
                end

                p_ttests = nan(1, numel(variablesToPlot));
                if ~isnan(p_interaction) && p_interaction < 0.05
                    for r = 1:numel(variablesToPlot)
                        [~, pval] = ttest2(g1(:,r), g2(:,r));
                        p_ttests(r) = pval;
                    end
                end

                metricOut = metric;
                if strcmp(metricOut,'Accuracy')
                    metricOut = 'Jaccard';
                end

                % Report F(df1, df2) and p for each within-subject effect
                fprintf('\n%s: repeated-measures ANOVA (LH/RH averaged, %s %s)\n', metricOut, cohort, ageGroup);
                fprintf('  Group:       F(%g, %g) = %.3f, p = %.4g\n', df1_group, df2_group, F_group, p_group);
                fprintf('  ROI:         F(%g, %g) = %.3f, p = %.4g\n', df1_roi,   df2_roi,   F_roi,   p_roi);
                fprintf('  Group x ROI: F(%g, %g) = %.3f, p = %.4g\n', df1_inter, df2_inter, F_inter, p_interaction);

                resultsSummary(end+1, :) = { ...
                    'both', metricOut, ...
                    p_group,       F_group, df1_group, df2_group, ...
                    p_roi,         F_roi,   df1_roi,   df2_roi, ...
                    p_interaction, F_inter, df1_inter, df2_inter, ...
                    p_ttests}; %#ok<AGROW>
            catch ME
                warning('Stats failed for %s (hemi-averaged): %s', metric, ME.message);
            end
        end
    end

    if strcmp(cohort,'ReplicationCombined')
        cohortLabel = 'Replication (adults + children)';
    else
        cohortLabel = sprintf('%s %s', cohort, ageGroup);
    end

    ax1.FontSize  = 26;
    ax1.LineWidth = 3;

    xticks(ax1, 1:numel(variablesToPlot));
    xticklabels(ax1, variablesToPlot);

    xlabel(ax1, 'ROI',   'FontSize', 20, 'FontWeight', 'bold');
    ylabel(ax1, 'Score', 'FontSize', 20, 'FontWeight', 'bold');

    legend(ax1, legendHandles, legendLabels, 'Location', 'best', 'FontSize', 17);
    grid(ax1,'off');

    % ================= SUBPLOT 2: double violins (spread across regions) =================
    % One double-violin per metric: left half = group1 (anat), right half =
    % group2 (prob). Each subject contributes one value = mean of the metric
    % across violinRegions (V3ab..IPS3), hemisphere-averaged. Colour encodes
    % the metric (matching subplot 1); side + shading encode anat vs prob.
    ax2 = axes('Position', axPosR);
    hold(ax2,'on');

    violinMetrics = violinOrder{plotIdx};

    for v = 1:numel(violinMetrics)
        metricV = violinMetrics{v};
        xpos    = violinXPositions(v);

        if ~isKey(metricHex, metricV)
            error('No color defined for metric "%s". Add it to metricHex.', metricV);
        end
        vColor = hex2rgb(metricHex(metricV));

        % hemi-averaged (subj x nRegions), then average across regions -> subj x 1
        Ma = loadMetricHemiAveraged(metricV, group1, violinRegions, false);
        Mp = loadMetricHemiAveraged(metricV, group2, violinRegions, false);
        va = mean(Ma, 2, 'omitnan');   % anat  (group1)
        vp = mean(Mp, 2, 'omitnan');   % prob  (group2)

        localHalfViolin(ax2, xpos, va, 'left',  vColor, 0.55);   % anat on the left
        localHalfViolin(ax2, xpos, vp, 'right', vColor, 0.25);   % prob on the right
    end

    % Merge with subplot 1: this panel carries NO axis lines, ticks, labels,
    % or title of its own. Convention: left half = anat (group1), right half
    % = prob (group2). x-positions are only used to space the three violins.
    box(ax2,'off');
    grid(ax2,'off');
    ax2.XAxis.Visible = 'off';        % no x ruler / ticks / labels
    ax2.YAxis.Visible = 'off';        % no y ruler / ticks / labels (shares ax1's)
    ax2.Color = 'none';               % transparent -> continuous with ax1
    xlim(ax2, [violinXPositions(1)-0.5, violinXPositions(end)+0.5]);

    % ---- shared y-axis across both subplots ----
    linkaxes([ax1 ax2], 'y');

    % Optional anat/prob key (delete these 5 lines for a fully bare panel):
    hAnat = patch(ax2, NaN, NaN, [0.35 0.35 0.35], 'FaceAlpha', 0.55, 'EdgeColor','none');
    hProb = patch(ax2, NaN, NaN, [0.35 0.35 0.35], 'FaceAlpha', 0.25, 'EdgeColor','none');
    legend(ax2, [hAnat hProb], ...
        {sprintf('%s (left)', prettyG1), sprintf('%s (right)', prettyG2)}, ...
        'Location','northeast', 'FontSize', 13, 'Box','off');

    titlesForPlot = metricLabels{plotIdx};
    titlesForPlot = strrep(strjoin(titlesForPlot, ', '), 'Accuracy', 'Jaccard');

    %sgtitle(sprintf('LH/RH Averaged Comparison (%s): %s', cohortLabel, titlesForPlot), ...
    %    'FontSize', 20, 'FontWeight', 'bold');
    % ---- VECTOR EXPORT ----
    fig = gcf;
    fig.Units = 'inches';
    fig.Position = [1 1 12 6];   % 2:1 to match the merged on-screen layout

    outName = sprintf('Compare_%s_%s_plot%d.pdf', cohort, ageGroup, plotIdx);
    exportgraphics(fig, outName, 'ContentType', 'vector');
end

%% === FORMAT OUTPUT (NO FILE WRITING) ===
if ~runStats || isempty(resultsSummary)
    resultsTable = [];
    return;
end

nRows      = size(resultsSummary, 1);
ttestCol   = 15;   % variable-length per-ROI t-test vector lives in this column
maxPostHoc = max(cellfun(@(x) numel(x), resultsSummary(:,ttestCol)));

for i = 1:nRows
    row = resultsSummary{i,ttestCol};
    resultsSummary{i,ttestCol} = [row, nan(1, maxPostHoc - numel(row))];
end

baseCols    = ttestCol - 1;   % 14 fixed columns before the per-ROI t-tests
posthocCols = maxPostHoc;

resultsExpanded = cell(nRows, baseCols + posthocCols);
resultsExpanded(:, 1:baseCols) = resultsSummary(:, 1:baseCols);

for i = 1:nRows
    rowVec = resultsSummary{i,ttestCol};
    for j = 1:posthocCols
        resultsExpanded{i, baseCols + j} = rowVec(j);
    end
end

posthocHeaders = arrayfun(@(i) sprintf('p_ROI%d_ttest', i), ...
    1:posthocCols, 'UniformOutput', false);

resultsTable = cell2table(resultsExpanded, 'VariableNames', ...
    [{'Hemisphere','Metric', ...
      'p_Group','F_Group','df1_Group','df2_Group', ...
      'p_ROI','F_ROI','df1_ROI','df2_ROI', ...
      'p_Group_x_ROI','F_Group_x_ROI','df1_Group_x_ROI','df2_Group_x_ROI'}, ...
     posthocHeaders])

%% === Helper: one half-violin (kernel density) + jittered points + mean ===
    function localHalfViolin(ax, xCenter, vals, side, faceColor, faceAlpha)
        % Draws a half violin on 'ax' at x = xCenter.
        %   side      : 'left' or 'right'
        %   faceColor : RGB (0..1); faceAlpha : patch transparency
        % Self-contained (uses ksdensity) so it does not depend on the
        % third-party violinplot version/API.
        maxWidth = 0.35;                 % max horizontal extent of the half violin
        if strcmpi(side,'left'), s = -1; else, s = 1; end

        vals = vals(:);
        vals = vals(~isnan(vals));
        if isempty(vals), return; end

        % --- density patch (needs >=2 distinct values) ---
        if numel(vals) >= 2 && (max(vals) > min(vals))
            [f, xi] = ksdensity(vals);
            f  = f / max(f) * maxWidth;                 % scale width
            px = xCenter + s * f(:);
            polyX = [px;      repmat(xCenter, numel(xi), 1)];
            polyY = [xi(:);   flipud(xi(:))];
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

end

