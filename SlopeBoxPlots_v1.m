%% =========================
% Prepare data for box plot
% =========================
load('Profile_Regression_Results_all_profiles_restrictive.mat')

Regression = string({Results.Regression});
Accepted   = string({Results.accepted_model});
Sector     = string({Results.Sector});
Slope      = [Results.Deming_slope];

keep = Regression == "NO3_vs_SA" & ...
       Accepted   == "linear"    & ...
       Slope ~= 0               & ...
       contains(Sector,"8090");

SlopeSel  = Slope(keep);
SectorSel = Sector(keep);

% Reduce sector names (remove latitude band)
BaseSector = erase(SectorSel, ["7080N","8090N"]);
BaseSector = strtrim(BaseSector);
BaseSector = categorical(BaseSector);

%% =========================
% Box plot of slopes by sector
% =========================

figure('Color','w','Position',[100 100 900 450]);

boxplot(SlopeSel, BaseSector, ...
    'Symbol','', ...
    'Whisker',1.5);

ylabel('Deming slope (Nitrate versus SA)')
title('High-latitude (80–90°N) Nitrate versus SA slopes by sector')
grid on
box on

%% =========================
% Overlay sector means
% =========================

hold on

Sectors = categories(categorical(BaseSector));
nsec = numel(Sectors);

meanSlopes = nan(nsec,1);

for i = 1:nsec
    meanSlopes(i) = mean(SlopeSel(BaseSector == Sectors{i}), 'omitnan');
end

plot(1:nsec, meanSlopes, 'kd', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',6);

%legend({'Mean'}, 'Location','best')

%% =========================
% Overlay 5–95 percentiles
% =========================

%p5  = nan(nsec,1);
%p95 = nan(nsec,1);

%for i = 1:nsec
%    vals = SlopeSel(BaseSector == Sectors{i});
%    p5(i)  = prctile(vals,5);
 %   p95(i) = prctile(vals,95);
%end

%for i = 1:nsec
%    plot([i i], [p5(i) p95(i)], 'r-', 'LineWidth',1.5)
%end

%legend({'Outliers','Mean',}, 'Location','best')