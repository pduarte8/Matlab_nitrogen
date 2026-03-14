%% =========================
% Profile-level regression analysis (Ice Regime)
% Adaptive linear diagnostics for small-N profiles
% =========================

clearvars; close all;

%% ------------------------
% Configuration
%% ------------------------
MinN = 5;           % minimum points for any regression
MinN_quad = 5;      % minimum points for quadratic regression
UpperDepthLimit = 100;
lambda = 1;         % Deming error ratio
Nboot = 200;        % bootstrap samples for quadratic diagnostics

%% ------------------------
% Load static lat/lon grids (once)
%% ------------------------
latGrid = single(geotiffread('O:\Norwegian Polar Institute\Papers\Nitrogen\ice_coverage_masks_Arctic_1980-2025\lat_5km.tif'));
lonGrid = single(geotiffread('O:\Norwegian Polar Institute\Papers\Nitrogen\ice_coverage_masks_Arctic_1980-2025\lon_5km.tif'));
[rows, cols] = size(latGrid);
gridIndex = reshape(1:rows*cols, rows, cols);
latVec = double(latGrid(:));
lonVec = double(lonGrid(:));
gridVec = double(gridIndex(:));
F_index = scatteredInterpolant(lonVec, latVec, gridVec, 'nearest');

%% ------------------------
% Data sets
%% ------------------------
DataSets = {
    'Barents_and_Arctic_European_Sector.mat', {'BarentsSector','EuropeanSector'}; 
    'Siberian_Sector.mat',                   {'SiberianSector7080N','SiberianSector8090N'}; 
    'Pacific_Sector.mat',                    {'PacificSector7080N','PacificSector8090N'}; 
    'Canadian_Sector.mat',                   {'CanadianSector7080N','CanadianSector8090N'}
};

%% ------------------------
% Preallocate Results
%% ------------------------
TotalProfiles = 20000;  % safe overestimate

Results(TotalProfiles) = struct( ...
    'ProfileID', NaN, ...
    'Sector', "", ...
    'Regression', "", ...
    'N', NaN, ...
    'AICc_linear', NaN, ...
    'AICc_quad', NaN, ...
    'quad_preferred', false, ...
    'accept_quadratic_AIC_only', false, ...
    'accepted_model', "none", ...
    'Deming_slope', NaN, ...
    'Deming_intercept', NaN, ...
    'quad_intercept', NaN, ...
    'quad_linear', NaN, ...
    'quad_curvature', NaN, ...
    'Lin', struct(), ...
    'Quad', struct(), ...
    'IceRegime', NaN ...
);

jj = 1;  % global counter

%% ------------------------
% Loop over data sets
%% ------------------------
for d = 1:size(DataSets,1)

    FileName = DataSets{d,1};
    Sectors  = DataSets{d,2};
    load(FileName);   % loads struct S

    for s = 1:numel(Sectors)

        Sector = string(Sectors{s});

        lonProf = S.(Sector).longitude.Value;
        latProf = S.(Sector).latitude.Value;
        UNIXSECS = S.(Sector).UNIXSECS.Value;
        p       = S.(Sector).PRESSURE.Value;
        depth   = S.(Sector).DEPTH.Value;
        SALT    = gsw_SA_from_SP(S.(Sector).SALNTY.Value,p,lonProf,latProf);
        NO3     = S.(Sector).NITRAT.Value;
        AOU     = S.(Sector).AOU.Value;
        WM      = S.(Sector).watermass.Value;
        [yearProf Month Day] = unixsecs2date_from_any_initial_year(UNIXSECS,1970); 
        
        [~, nprof] = size(SALT);

        for j = 1:nprof

            % --- Atlantic Water selection
            k = find(WM(:,j) >= 30 & WM(:,j) <= 50 & depth(:,j) >= UpperDepthLimit);
            if numel(k) < MinN
                continue
            end

            % --- Assign ice regime

            yr = yearProf(j);
            if yr >= 1980
                maskFile = sprintf('O:/Norwegian Polar Institute/Papers/Nitrogen/ice_coverage_masks_Arctic_1980-2025/seaice_frequency_max_%d_5km.tif', yr);
                maskGrid = readgeoraster(maskFile);
                maskGrid = logical(maskGrid);
                idxMask = F_index(lonProf(j), latProf(j));
                IceRegime = maskGrid(idxMask);  % 1 = persistent ice, 0 = open water
    
                %% ======================================================
                % NO3 vs Salinity (same as before)
                %% ======================================================
                x = SALT(k,j);
                y = NO3(k,j);
                idx = isfinite(x) & isfinite(y);
                x = x(idx); y = y(idx);
                N = numel(x);
                if N < MinN
                    continue
                end
                Lin = diagnose_linearity_adaptive_restrictive(x,y);
                if N >= MinN_quad
                    Quad = diagnose_quadratic_structure(x,y);
                else
                    Quad = struct('passes_complementary_tests',false,'quad',struct());
                end
                [b1,b0] = deming_regress(x,y,lambda);
                yhatL = b0 + b1*x; RSS_L = sum((y - yhatL).^2); kL=2;
                AICc_L = N*log(RSS_L/N) + 2*kL + (2*kL*(kL+1))/(N-kL-1);
                if N >= MinN_quad
                    xc = x - mean(x);
                    Xq = [ones(N,1) xc xc.^2];
                    betaQ = Xq \ y;
                    yhatQ = Xq * betaQ; RSS_Q = sum((y - yhatQ).^2); kQ=3;
                    AICc_Q = N*log(RSS_Q/N) + 2*kQ + (2*kQ*(kQ+1))/(N-kQ-1);
                    quad_preferred = (AICc_Q - AICc_L) < -2;
                    accept_quadratic_all = quad_preferred && Quad.passes_complementary_tests;
                    accept_quadratic_AIC_only = quad_preferred && ~Quad.passes_complementary_tests;
                else
                    betaQ = [NaN; NaN; NaN]; AICc_Q=NaN; quad_preferred=false;
                    accept_quadratic_all=false; accept_quadratic_AIC_only=false;
                end
                if accept_quadratic_all
                    accepted_model = "quadratic";
                elseif Lin.is_linear
                    accepted_model = "linear";
                elseif accept_quadratic_AIC_only
                    accepted_model = "quadratic";
                else
                    accepted_model = "none";
                end
                % --- Store results
                Results(jj) = struct( ...
                    'ProfileID', j, ...
                    'Sector', Sector, ...
                    'Regression', 'NO3_vs_SA', ...
                    'N', N, ...
                    'AICc_linear', AICc_L, ...
                    'AICc_quad', AICc_Q, ...
                    'quad_preferred', quad_preferred, ...
                    'accept_quadratic_AIC_only', accept_quadratic_AIC_only, ...
                    'accepted_model', accepted_model, ...
                    'Deming_slope', b1, ...
                    'Deming_intercept', b0, ...
                    'quad_intercept', betaQ(1), ...
                    'quad_linear', betaQ(2), ...
                    'quad_curvature', betaQ(3), ...
                    'Lin', Lin, ...
                    'Quad', Quad, ...
                    'IceRegime', IceRegime ...
                );
                jj = jj + 1;
    
                %% ======================================================
                % AOU vs NO3 (same structure, add IceRegime)
                %% ======================================================
                x = NO3(k,j);
                y = AOU(k,j);
                idx = isfinite(x) & isfinite(y);
                x = x(idx); y = y(idx);
                N = numel(x);
                if N < MinN
                    continue
                end
                Lin = diagnose_linearity_adaptive_restrictive(x,y);
                if N >= MinN_quad
                    Quad = diagnose_quadratic_structure(x,y);
                else
                    Quad = struct('passes_complementary_tests',false,'quad',struct());
                end
                [b1,b0] = deming_regress(x,y,lambda);
                yhatL = b0 + b1*x; RSS_L = sum((y - yhatL).^2);
                AICc_L = N*log(RSS_L/N) + 2*kL + (2*kL*(kL+1))/(N-kL-1);
                if N >= MinN_quad
                    xc = x - mean(x); Xq = [ones(N,1) xc xc.^2]; betaQ = Xq \ y;
                    yhatQ = Xq * betaQ; RSS_Q = sum((y - yhatQ).^2);
                    AICc_Q = N*log(RSS_Q/N) + 2*kQ + (2*kQ*(kQ+1))/(N-kQ-1);
                    quad_preferred = (AICc_Q - AICc_L) < -2;
                    accept_quadratic_all = quad_preferred && Quad.passes_complementary_tests;
                    accept_quadratic_AIC_only = quad_preferred && ~Quad.passes_complementary_tests;
                else
                    betaQ = [NaN; NaN; NaN]; AICc_Q=NaN; quad_preferred=false;
                    accept_quadratic_all=false; accept_quadratic_AIC_only=false;
                end
                if accept_quadratic_all
                    accepted_model = "quadratic";
                elseif Lin.is_linear
                    accepted_model = "linear";
                elseif accept_quadratic_AIC_only
                    accepted_model = "quadratic";
                else
                    accepted_model = "none";
                end
                Results(jj) = struct( ...
                    'ProfileID', j, ...
                    'Sector', Sector, ...
                    'Regression', 'AOU_vs_NO3', ...
                    'N', N, ...
                    'AICc_linear', AICc_L, ...
                    'AICc_quad', AICc_Q, ...
                    'quad_preferred', quad_preferred, ...
                    'accept_quadratic_AIC_only', accept_quadratic_AIC_only, ...
                    'accepted_model', accepted_model, ...
                    'Deming_slope', b1, ...
                    'Deming_intercept', b0, ...
                    'quad_intercept', betaQ(1), ...
                    'quad_linear', betaQ(2), ...
                    'quad_curvature', betaQ(3), ...
                    'Lin', Lin, ...
                    'Quad', Quad, ...
                    'IceRegime', IceRegime ...
                );
                jj = jj + 1;
            end
        end
    end
end

%% ------------------------
% Trim and save
%% ------------------------
Results = Results(1:jj-1);
save('Profile_Regression_Results_all_profiles_ice_regime.mat','Results','-v7.3');