%% ============================================================
% Frequency distributions of linear slopes (profile-level)
% Uses Results from Nitrate_vs_Salinity_vs_AOU_patterns_profile_data_v4
% ============================================================

clearvars; close all;

%% ------------------------
% Load results
%% ------------------------
load('Profile_Regression_Results_all_profiles_restrictive.mat','Results');

%% ------------------------
% Configuration
%% ------------------------
FigDir = 'Figures_Profile_Slope_Distributions';
if ~exist(FigDir,'dir')
    mkdir(FigDir)
end

VarPairs = {'NO3_vs_SA','AOU_vs_NO3'};

%% ------------------------
% Extract metadata helpers
%% ------------------------
Sectors = unique(string({Results.Sector}));

get_latband = @(s) infer_latband(s);



%% ------------------------
% Loop over variable pairs, sectors, latitude bands
%% ------------------------
for v = 1:numel(VarPairs)
    regname = VarPairs{v};

    for s = 1:numel(Sectors)

        sector = Sectors(s);
        latband = get_latband(sector);

        if latband == ""
            continue
        end

        %% ------------------------
        % Select valid linear profiles
        %% ------------------------
        idx = ...
            string({Results.Regression}) == regname & ...
            string({Results.Sector})     == sector & ...
            string({Results.accepted_model}) == "linear" & ...
            isfinite([Results.Deming_slope]);

        slopes = [Results(idx).Deming_slope];

        if numel(slopes) < 10
            continue
        end

        %% ------------------------
        % Statistics
        %% ------------------------
        med_slope = median(slopes);
        mean_slope = mean(slopes);

        %% ------------------------
        % Plot
        %% ------------------------
        fig = figure('Visible','on'); hold on;

        histogram(slopes, ...
            'Normalization','count', ...
            'NumBins', round(sqrt(numel(slopes))), ...
            'FaceColor',[0.3 0.3 0.8], ...
            'EdgeColor','none');

        yl = ylim;

        plot([med_slope med_slope], yl, ...
            'r','LineWidth',2, ...
            'DisplayName','Median');

        plot([mean_slope mean_slope], yl, ...
            '--k','LineWidth',1.5, ...
            'DisplayName','Mean');

        %% ------------------------
        % Labels & title
        %% ------------------------
        xlabel('Deming slope','FontSize',11)
        ylabel('Frequency','FontSize',11)

        title(sprintf('%s | %s | %s\nN = %d profiles', ...
            regname, sector, latband, numel(slopes)), ...
            'Interpreter','none','FontSize',11)

        legend('Location','best')
        grid on
        box on

        %% ------------------------
        % Save
        %% ------------------------
        fname = sprintf('%s_%s_%s',regname,sector,latband);
        fname = strrep(fname,'–','_');
        fname = strrep(fname,' ','_');

        saveas(fig, fullfile(FigDir,[fname '.png']));
        savefig(fig, fullfile(FigDir,[fname '.fig']));
        close(fig)

    end
end

function latband = infer_latband(sector)
    if contains(sector,'7080')
        latband = "70–80N";
    elseif contains(sector,'8090')
        latband = "80–90N";
    else
        latband = "";
    end
end