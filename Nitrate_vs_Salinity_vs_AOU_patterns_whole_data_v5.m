clear all

%% =========================
% INPUT
%% =========================

FigDir = 'Figures_whole_sector_restrictive';
if ~exist(FigDir,'dir')
    mkdir(FigDir)
end

%% =========================
% Axis labels
%% =========================
LabelMap = struct();

LabelMap.NO3_vs_SA.x = 'Salinity absolute';
LabelMap.NO3_vs_SA.y = 'Nitrate (\mumol kg^{-1})';

LabelMap.AOU_vs_SA.x = 'Salinity absolute';
LabelMap.AOU_vs_SA.y = 'AOU (\mumol kg^{-1})';

LabelMap.AOU_vs_NO3.x = 'Nitrate (\mumol kg^{-1})';
LabelMap.AOU_vs_NO3.y = 'AOU (\mumol kg^{-1})';

DataSets = {
    'Barents_and_Arctic_European_Sector.mat', {'BarentsSector','EuropeanSector'};
    'Siberian_Sector.mat',                   {'SiberianSector7080N','SiberianSector8090N'};
    'Pacific_Sector.mat',                    {'PacificSector7080N','PacificSector8090N'};
    'Canadian_Sector.mat',                   {'CanadianSector7080N','CanadianSector8090N'};
};

j = 1;
Results = struct();
UpperDepthLimit = 100;

for d = 1:size(DataSets,1)

    FileName = DataSets{d,1};
    Sectors  = DataSets{d,2};
    load(FileName);

    for s = 1:length(Sectors)

        Sector = string(Sectors{s});

        %% -------------------------
        % Load sector data
        %% -------------------------
        depth = S.(Sector).DEPTH.Value;
        lon   = S.(Sector).longitude.Value;
        lat   = S.(Sector).latitude.Value;
        SP    = S.(Sector).SALNTY.Value;
        p     = S.(Sector).PRESSURE.Value;
        SA    = gsw_SA_from_SP(SP,p,lon,lat);
        NO3   = S.(Sector).NITRAT.Value;
        AOU   = S.(Sector).AOU.Value;
        WM    = S.(Sector).watermass.Value;

        %% Atlantic Water selection
        [ii,jj] = find(WM >= 30 & WM <= 50 & depth >= UpperDepthLimit);

        SA_v  = SA(sub2ind(size(SA),ii,jj));
        NO3_v = NO3(sub2ind(size(NO3),ii,jj));
        AOU_v = AOU(sub2ind(size(AOU),ii,jj));

        %% =========================
        % Regression pairs
        %% =========================
        Pairs = {
            SA_v,  NO3_v, 'NO3_vs_SA';
            SA_v,  AOU_v, 'AOU_vs_SA';
            NO3_v, AOU_v, 'AOU_vs_NO3'
        };

        for r = 1:size(Pairs,1)

            x = Pairs{r,1};
            y = Pairs{r,2};
            reg_name = Pairs{r,3};

            %% -------------------------
            % Clean
            %% -------------------------
            idx = isfinite(x) & isfinite(y);
            x = x(idx);
            y = y(idx);
            N = numel(x);
            if N < 6
                continue
            end

            %% =========================
            % Linear diagnostics
            %% =========================
            Lin = diagnose_linearity_adaptive_restrictive(x,y);

            %% =========================
            % Quadratic diagnostics
            %% =========================
            Quad = diagnose_quadratic_structure(x,y); % only fit + metrics

            %% =========================
            % Compute AIC for linear and quadratic fits
            %% -------------------------
            % Linear Deming
            [b1,b0] = deming_regress(x,y);
            yhatL = b0 + b1*x;
            RSS_L = sum((y - yhatL).^2);
            kL = 2;
            AICc_L = N*log(RSS_L/N) + 2*kL + (2*kL*(kL+1))/(N-kL-1);

            % Quadratic OLS (centered)
            xc = x - mean(x);
            Xq = [ones(N,1) xc xc.^2];
            betaQ = Xq\y;
            yhatQ = Xq*betaQ;
            RSS_Q = sum((y - yhatQ).^2);
            kQ = 3;
            AICc_Q = N*log(RSS_Q/N) + 2*kQ + (2*kQ*(kQ+1))/(N-kQ-1);

            %% -------------------------
            % --- Determine quadratic acceptance
            accept_quadratic_all = (AICc_Q < AICc_L - 2) && Quad.passes_complementary_tests;
            accept_quadratic_AIC_only = (AICc_Q < AICc_L - 2) && ~Quad.passes_complementary_tests;
            
            % --- Determine final accepted model
            if accept_quadratic_all
                accepted_model = "quadratic";
            elseif accept_quadratic_AIC_only
                accepted_model = "quadratic (AIC only)";
            elseif Lin.is_linear
                accepted_model = "linear";
            else
                accepted_model = "none";
            end
            
            % --- Store in Results
            Results(j).Sector                     = Sector;
            Results(j).Regression                 = reg_name;
            Results(j).N                          = N;
            Results(j).AICc_linear                = AICc_L;
            Results(j).AICc_quad                  = AICc_Q;
            Results(j).accept_quadratic_all       = accept_quadratic_all;
            Results(j).accept_quadratic_AIC_only  = accept_quadratic_AIC_only;
            Results(j).accepted_model             = accepted_model;
            
            % ---- Linear diagnostics
            Results(j).Lin                        = Lin;
            Results(j).Deming_slope               = b1;
            Results(j).Deming_intercept           = b0;
            
            % ---- Quadratic diagnostics
            Results(j).Quad                       = Quad;
            Results(j).quad_intercept             = Quad.quad.intercept;
            Results(j).quad_linear                = Quad.quad.linear;
            Results(j).quad_curvature             = Quad.quad.curvature;

            Results(j).File = FileName;

            %% =========================
            % FIGURE: scatter + accepted regression
            %% =========================
            fig = figure('Visible','on'); hold on

            scatter(x,y,6,'k','filled','MarkerFaceAlpha',0.25)

            % ---- Linear fit
            if Lin.is_linear
                xx = linspace(min(x),max(x),200);
                yy = b0 + b1*xx;
                plot(xx,yy,'b','LineWidth',2,'DisplayName','Linear');
            end

            % ---- Quadratic fit
            xx = linspace(min(x), max(x), 200);
            xx0 = xx - mean(x);  % center for plotting
            
            if accept_quadratic_all
                yy = betaQ(1) + betaQ(2)*xx0 + betaQ(3)*xx0.^2;
                plot(xx, yy, 'r', 'LineWidth', 2, 'DisplayName', 'Quadratic');
            elseif accept_quadratic_AIC_only
                yy = betaQ(1) + betaQ(2)*xx0 + betaQ(3)*xx0.^2;
                plot(xx, yy, 'r--', 'LineWidth', 2, 'DisplayName', 'Quadratic (AIC only)');
            end

            xlabel(LabelMap.(reg_name).x,'FontSize',11)
            ylabel(LabelMap.(reg_name).y,'FontSize',11)
            title(Sector,'Interpreter','none','FontSize',12)

            grid on
            box on
            legend('Location','best')

            % ---- Save figure
            sector_str = char(Sector);
            file_str = fullfile(FigDir, [sector_str, '_', reg_name]);
            savefig(fig, [file_str, '.fig']);
            saveas(fig, [file_str, '.png']);
            close(fig)

            j = j + 1;

        end
    end
end