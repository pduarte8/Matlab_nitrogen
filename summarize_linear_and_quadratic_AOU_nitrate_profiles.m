clearvars;

%load('Profile_Regression_Results_all_profiles_v4.mat','Results');
load('Profile_Regression_Results_all_profiles_restrictive.mat','Results');
%% =========================
% Select NO3 vs SA regressions
% =========================
Regression = string({Results.Regression});
Accepted   = string({Results.accepted_model});
Sector     = string({Results.Sector});


is_AOU_NO3 = Regression == "AOU_vs_NO3";

R = Results(is_AOU_NO3);

Accepted = Accepted(is_AOU_NO3);
Sector   = Sector(is_AOU_NO3);

%% =========================
% Latitude band from sector name
% =========================
LatBand = strings(size(Sector));
LatBand(contains(Sector,"7080")) = "70–80N";
LatBand(contains(Sector,"8090")) = "80–90N";
LatBand(LatBand=="") = "Other";

GroupLabel = extractBetween(Sector, 1, strlength(Sector) - 5) + " " + LatBand;
Groups = unique(GroupLabel,'stable');
newOrder = [1,3,5,7,2,4,6,8];
Groups = Groups(newOrder);
%% =========================
% Initialize counters
% =========================
nGroups = numel(Groups);

n_lin_pos   = zeros(nGroups,1);
n_lin_neg   = zeros(nGroups,1);
n_quad_src  = zeros(nGroups,1);
n_quad_sink = zeros(nGroups,1);
n_none      = zeros(nGroups,1);

%% =========================
% Loop over groups
% =========================
for g = 1:nGroups
    ig = GroupLabel == Groups(g);

    % ----- Linear
    il = Accepted(ig) == "linear";
    if any(il)
        slopes = [R(ig).Deming_slope];
        slopes = slopes(il);

        n_lin_pos(g) = sum(slopes > 0);
        n_lin_neg(g) = sum(slopes < 0);
    end

    % ----- Quadratic
    iq = Accepted(ig) == "quadratic";
    if any(iq)
        b = [R(ig).quad_linear];
        c = [R(ig).quad_curvature];

        b = b(iq);
        c = c(iq);

        n_quad_src(g)  = sum(b > 0 & c < 0);
        n_quad_sink(g) = sum(b < 0 & c > 0);
        n_quad_other(g) = sum(b > 0 & c > 0) + sum(b < 0 & c < 0);
    end

    % ----- None
    n_none(g) = sum(Accepted(ig) == "none") + n_quad_other(g);
end

%% =========================
% Plot
% =========================
figure; hold on
GroupCategories = categorical(Groups, Groups, 'Ordinal', true);
Total = n_lin_pos+n_lin_neg+n_quad_src+n_quad_sink+n_none;
bar(GroupCategories, ...
    [n_lin_pos n_lin_neg n_quad_src n_quad_sink n_none], ...
    'stacked');

ylabel('Number of profiles');

legend({ ...
    'Linear (positive slope)', ...
    'Linear (negative slope)', ...
    'Quadratic (source-like)', ...
    'Quadratic (sink-like)', ...
    'None'}, ...
    'Location','best');
    title('Profile-level AOU versus nitrate regression outcomes');

grid on
box on
xtickangle(35);