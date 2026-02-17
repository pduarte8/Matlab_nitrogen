clearvars;

load('Profile_Regression_Results_all_profiles_v4.mat','Results');

%% =========================
% Selection
% =========================

% --- NO3 vs SA only
is_NO3_SA = strcmp({Results.Regression}, 'NO3_vs_SA');

% --- Accepted quadratic models
accepted_model = string({Results.accepted_model});
is_quad = (accepted_model == "quadratic");

% --- OPTIONAL: require complementary tests
require_complementary = false;

if require_complementary
    pass_comp = arrayfun(@(r) ...
        isfield(r.Quad,'passes_complementary_tests') && r.Quad.passes_complementary_tests, ...
        Results);
else
    pass_comp = true(size(Results));
end

keep = is_NO3_SA & is_quad; %& pass_comp;

R = Results(keep);

fprintf('Total significant quadratic NO3–SA profiles: %d\n', numel(R));

%% =========================
% Classify quadratic types
% =========================

b = [R.quad_linear];
c = [R.quad_curvature];

is_source_like = (b > 0) & (c < 0);  % source
is_sink_like   = (b < 0) & (c > 0);  % sink

% Sanity check
unclassified = ~(is_source_like | is_sink_like);

fprintf('Source-like quadratics: %d\n', sum(is_source_like));
fprintf('Sink-like quadratics:   %d\n', sum(is_sink_like));
fprintf('Other / mixed signs:    %d\n', sum(unclassified));

%% =========================
% Sector and latitude parsing
% =========================

sectors = string({R.Sector});

% Latitude band inference from sector name
latband = strings(size(sectors));
latband(contains(sectors,'7080')) = "70–80N";
latband(contains(sectors,'8090')) = "80–90N";
latband(latband=="") = "Other";

%% =========================
% Tabulate frequencies
% =========================

T = table( ...
    sectors(:), latband(:), ...
    is_source_like(:), is_sink_like(:), ...
    'VariableNames', {'Sector','LatBand','SourceLike','SinkLike'});

Summary = groupsummary(T, {'Sector','LatBand'}, 'sum', {'SourceLike','SinkLike'});

Summary.Total = Summary.sum_SourceLike + Summary.sum_SinkLike;

Summary.Frac_SourceLike = Summary.sum_SourceLike ./ Summary.Total;
Summary.Frac_SinkLike   = Summary.sum_SinkLike   ./ Summary.Total;

disp(Summary);

%% =========================
% Visualization
% =========================

figure; hold on

cats = categorical(Summary.Sector + " (" + Summary.LatBand + ")");
bar(cats, [Summary.sum_SourceLike Summary.sum_SinkLike], 'stacked');

ylabel('Number of profiles');
legend({'Source-like (b>0, c<0)','Sink-like (b<0, c>0)'}, 'Location','best');
title('Quadratic NO_3–SA regressions (profiles)');
grid on
box on