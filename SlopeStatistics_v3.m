%% =========================================================
% NO3–SA slope statistics & sector comparison (80–90N)
% =========================================================

clearvars; clc;
alpha = 0.05;
load('Profile_Regression_Results_all_profiles_v4.mat','Results');

%% =========================
% Selection
% =========================

Regression = string({Results.Regression});
Accepted   = string({Results.accepted_model});
SectorStr  = string({Results.Sector});
Slope      = [Results.Deming_slope];

keep = Regression == "NO3_vs_SA" & ...
       Accepted   == "linear"    & ...
       contains(SectorStr,"8090");

R = Results(keep);

Sector = string({R.Sector});
Slope  = [R.Deming_slope];

%% =========================
% Base sector names
% =========================

BaseSector = erase(Sector, ["7080N","8090N"]);
BaseSector = strtrim(BaseSector);
BaseSector = categorical(BaseSector);

Sectors = unique(BaseSector,'stable');
nSec = numel(Sectors);

%% =========================
% Descriptive statistics
% =========================

N       = zeros(nSec,1);
Mean    = zeros(nSec,1);
Median  = zeros(nSec,1);
P2p5    = zeros(nSec,1);
P97p5   = zeros(nSec,1);
CI_low  = zeros(nSec,1);
CI_high = zeros(nSec,1);

nboot = 5000;

for i = 1:nSec

    idx = BaseSector == Sectors(i);
    s   = Slope(idx);

    N(i)      = numel(s);
    Mean(i)   = mean(s,'omitnan');
    Median(i) = median(s,'omitnan');

    P2p5(i)  = prctile(s,2.5);
    P97p5(i) = prctile(s,97.5);

    % --- Bootstrap CI for the mean
    bootstat = bootstrp(nboot,@mean,s);
    CI = prctile(bootstat,[2.5 97.5]);

    CI_low(i)  = CI(1);
    CI_high(i) = CI(2);
end

%% =========================
% Assemble summary table
% =========================

SlopeStats = table( ...
    Sectors(:), N, ...
    Mean, Median, ...
    P2p5, P97p5, ...
    CI_low, CI_high, ...
    'VariableNames', { ...
        'Sector', 'N_profiles', ...
        'Mean_slope', 'Median_slope', ...
        'P2p5', 'P97p5', ...
        'CI95_low', 'CI95_high'});

disp(SlopeStats)

%% =========================
% Statistical tests
% =========================

% --- Kruskal–Wallis (global test)
[p_kw,~,stats_kw] = kruskalwallis(Slope, BaseSector, 'off');

fprintf('\nKruskal–Wallis test:\n');
fprintf('p = %.4g\n', p_kw);
disp("Group names used by Kruskal–Wallis:")
disp(stats_kw.gnames)
% --- Post-hoc Dunn–Šidák tests
if p_kw < 0.05
    comp = multcompare(stats_kw, ...
        'CType','dunn-sidak', ...
        'Display','off');
    pvals = comp(:,6);   % raw p-values
    [p_fdr] = mafdr(pvals, 'BHFDR', true);
    h_fdr = p_fdr < alpha;
    PostHoc = array2table(comp, ...
        'VariableNames',{'Group1','Group2','Lower','Diff','Upper','pValue'});

    GroupNames = string(stats_kw.gnames);

    PostHoc.Sector1 = GroupNames(PostHoc.Group1);
    PostHoc.Sector2 = GroupNames(PostHoc.Group2);

    PostHoc = movevars(PostHoc,{'Sector1','Sector2'},'Before','Lower');
    PostHoc.p_FDR = p_fdr;
    PostHoc.h_FDR = h_fdr;
    disp('Post-hoc Dunn–Šidák comparisons:')
    disp(PostHoc(:,{'Sector1','Sector2','pValue'}))
else
    disp('No significant global difference → post-hoc tests not applied.')
end