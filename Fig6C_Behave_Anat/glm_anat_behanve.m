% Plots VSB predictions from each anatomical measure 
% Each measure is one column: beta bars on top, predicted-vs-measured scatter below.
%
% Requires MATLAB 2025+

clear; clc;

%% ---- Setup ----
baseDataDir = '.';
funcDir     = fullfile(baseDataDir, 'func');
devTag      = 'adults';

measures    = {'ThickStd', 'MeanCurv', 'Depth'};
measLabels  = {'Thickness Std', 'Mean Curvature', 'Depth'};  % display names for titles

labelSource = 'func';
predMode    = 'meanlr';
labels      = {'V3ab','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5'};

%% ---- Load behavioral data ----
behFile = './observation_behavior.xlsx';
B       = readtable(behFile, "VariableNamingRule","preserve");

subBeh = string(B.Subject);
y_vsb  = zscore_omitnan(B.VSB);

%% ---- Figure layout (manual axes positioning) ----
nMeas = numel(measures);

% All sizes in inches
axW       = 2.0;    % axis width (same for both rows)
barH      = 5.0;    % tall rectangle for beta bars
scatH     = 2.0;    % smaller square for scatter
colGap    = 0.9;    % horizontal gap between columns
leftMarg  = 1.0;    % left margin (for y-axis labels)
rightMarg = 0.4;
topMarg   = 0.6;    % top margin
midGap    = 0.8;    % vertical gap between rows
botMarg   = 1.0;    % bottom margin (for x-axis labels)

figW = leftMarg + nMeas * axW + (nMeas-1) * colGap + rightMarg;
figH = topMarg + barH + midGap + scatH + botMarg;

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 figW figH]);

% Precompute axes positions [left bottom width height] in normalized units
axPos_bar  = cell(nMeas, 1);
axPos_scat = cell(nMeas, 1);

for mm = 1:nMeas
    left = (leftMarg + (mm-1) * (axW + colGap)) / figW;
    % bar: top row
    bot_bar  = (botMarg + scatH + midGap) / figH;
    axPos_bar{mm}  = [left, bot_bar,  axW/figW, barH/figH];
    % scatter: bottom row
    bot_scat = botMarg / figH;
    axPos_scat{mm} = [left, bot_scat, axW/figW, scatH/figH];
end

dataDir = funcDir;
prefix  = 'func';
cd(dataDir);

%% ---- Loop over measures ----
for mm = 1:nMeas

    whichMeasure = measures{mm};

    fL = sprintf('lh.%s.stats_table.%s.%s.mat', prefix, whichMeasure, devTag);
    fR = sprintf('rh.%s.stats_table.%s.%s.mat', prefix, whichMeasure, devTag);

    if ~exist(fL,'file') || ~exist(fR,'file')
        warning('Missing files for %s — skipping.', whichMeasure);
        continue;
    end

    SL = load(fL, 'dataTable');  SR = load(fR, 'dataTable');
    TL = SL.dataTable;           TR = SR.dataTable;

    subL = string(TL.subjects2);
    subR = string(TR.subjects2);

    subMaster = intersect(intersect(subL, subR, 'stable'), subBeh, 'stable');

    [~, idxBeh] = ismember(subMaster, subBeh);
    [~, idxL]   = ismember(subMaster, subL);
    [~, idxR]   = ismember(subMaster, subR);

    % Build predictor matrix
    X = nan(numel(subMaster), numel(labels));
    for ii = 1:numel(labels)
        region = labels{ii};
        if ~ismember(region, TL.Properties.VariableNames) || ...
           ~ismember(region, TR.Properties.VariableNames)
            warning('Region %s not found for %s. Filling NaNs.', region, whichMeasure);
            continue;
        end
        L = TL.(region)(idxL);
        R = TR.(region)(idxR);
        X(:,ii) = (L + R) ./ 2;
        X(~isfinite(X(:,ii)), ii) = NaN;
    end

    y_vsb_m = y_vsb(idxBeh);
    Xz      = zscore_cols_omitnan(X);

    % Fit LOO GLM
    [~, b_vsb, yfit_vsb] = fit_glm_nanrows_loo(Xz, y_vsb_m);

    % ---- Row 1: Beta bar chart ----
    ax1 = axes('Units', 'normalized', 'Position', axPos_bar{mm}); %#ok<LAXES>
    plot_beta_barh(b_vsb, labels);
    title(measLabels{mm}, 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');

    % ---- Row 2: Scatter ----
    ax2 = axes('Units', 'normalized', 'Position', axPos_scat{mm}); %#ok<LAXES>
    plot_pred_scatter(yfit_vsb, y_vsb_m, whichMeasure);

    xlabel(sprintf('Predicted\nVisuospatial Bias'), ...
        'FontSize', 12, 'FontWeight', 'bold');

    if mm == 1
        ylabel('\bfMeasured', 'FontSize', 12);
    else
        ylabel('');
    end

end

%% ---- Save PDF ----
pdfOut = fullfile(dataDir, sprintf('VSB_anatMeasures_%s_%s.pdf', labelSource, predMode));
exportgraphics(fig, pdfOut, 'ContentType', 'vector');
fprintf('Saved: %s\n', pdfOut);

%% ======================== Helper functions ========================

function z = zscore_omitnan(x)
    x  = x(:);
    mu = mean(x, 'omitnan');
    sg = std(x,  'omitnan');
    z  = (x - mu) ./ sg;
end

function Xz = zscore_cols_omitnan(X)
    Xz = X;
    for c = 1:size(X,2)
        Xz(:,c) = zscore_omitnan(X(:,c));
    end
end

function [good, b, yfit_full] = fit_glm_nanrows_loo(X, y)
    good = all(isfinite(X),2) & isfinite(y);
    Xg   = X(good,:);
    yg   = y(good);
    n    = numel(yg);

    yfit_g = nan(n,1);
    b_all  = nan(size(Xg,2)+1, n);

    for ii = 1:n
        trainIdx        = true(n,1);
        trainIdx(ii)    = false;
        b_fold          = glmfit(Xg(trainIdx,:), yg(trainIdx), 'normal');
        yfit_g(ii)      = glmval(b_fold, Xg(ii,:), 'identity');
        b_all(:,ii)     = b_fold;
    end

    yfit_full             = nan(size(y));
    yfit_full(find(good)) = yfit_g;
    b = mean(b_all, 2, 'omitnan');
end

function plot_beta_barh(b, labels)
    bh = barh(b(2:end));
    bh.FaceColor = 'flat';
    bh.CData     = b(2:end);

    colormap('nebula');
    caxis([-1.5 1.5]);
    xlim([-1.4 1.4]); %space of how far the plots exttend within the plots

    set(gca, ...
        'FontSize',    13, ...
        'TickDir',     'out', ...
        'Box',         'off', ...
        'YTickLabel',  labels, ...
        'LineWidth',   3,...
        'TickLength',  [0.02 0.02]);

    % thin gray zero line
    xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);

    xlabel('\beta', 'FontSize', 14, 'FontWeight', 'bold');
end

function plot_pred_scatter(yfit, y, labelStr)
    if nargin < 3; labelStr = ''; end
    good = isfinite(yfit) & isfinite(y);

    scatter(yfit(good), y(good), 20, 'filled', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k');
    hold on;
    l           = lsline;
    l.LineWidth = 2.5;
    l.Color     = 'k';
    hold off;

    % Symmetric axis limits, data-driven
    allVals = [yfit(good); y(good)];
    lim     = ceil(max(abs(allVals)) * 10) / 10 + 0.1;
    lim     = max(lim, 2.0);
    set(gca, ...
        'XLim',       [-lim lim], ...
        'YLim',       [-lim lim], ...
        'TickDir',    'out', ...
        'FontSize',   12, ...
        'Box',        'off', ...
        'LineWidth',   3,...
        'TickLength', [0.02 0.02]);

    % Correlation
    if nnz(good) >= 3
        [R, P] = corrcoef(yfit(good), y(good));
        r = R(1,2);  p = P(1,2);
    else
        r = NaN;  p = NaN;
    end

    fprintf('%s  VSB: N=%d, r=%.2f, p=%.3g\n', labelStr, nnz(good), r, p);

    % Annotation bottom-left
    ax = gca;
    x0 = ax.XLim(1) + 0.05 * range(ax.XLim);
    y0 = ax.YLim(1) + 0.15 * range(ax.YLim);

    if isfinite(r)
        txt = sprintf('r=%.2f, p=%.4g', r, p);
    else
        txt = sprintf('N=%d (insuff.)', nnz(good));
    end
    text(x0, y0, txt, 'FontSize', 10, 'Interpreter', 'none');
end
