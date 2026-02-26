%% =========================================================
% Representative high-latitude NO3–SA profiles (two-panel)
% Uses Deming regression + UNIX date
% =========================================================

clearvars; close all; clc;

%% =========================================================
% Load sector data (KEEP SEPARATE)
% =========================================================

load('Canadian_Sector.mat','S');      S_CAN = S;
load('Pacific_Sector.mat','S');       S_PAC = S;
load('Siberian_Sector.mat','S');      S_SIB = S;
load('Barents_and_Arctic_European_Sector.mat','S');
S_EUR = S;
clear S

%% =========================================================
% Load regression results
% =========================================================

load('Profile_Regression_Results_all_profiles_restrictive.mat','Results');

%% =========================================================
% Filter Results: NO3 vs SA, linear, positive slope, 80–90N
% =========================================================

Regression = string({Results.Regression});
Accepted   = string({Results.accepted_model});
SectorStr  = string({Results.Sector});
Slope      = [Results.Deming_slope];

keep = Regression == "NO3_vs_SA" & ...
       Accepted   == "linear"    & ...
       Slope > 0               & ...
       contains(SectorStr,"8090");

Rsel = Results(keep);

%% =========================================================
% Select one representative profile per BASE sector
% (max N_eff)
% =========================================================

SectorSel  = string({Rsel.Sector});
BaseSel    = erase(SectorSel, ["7080N","8090N"]);
UniqueBase = unique(BaseSel);

Selected = struct();

for i = 1:numel(UniqueBase)

    bs = UniqueBase(i);
    idx = BaseSel == bs;
    Rs  = Rsel(idx);

    if isempty(Rs)
        continue
    end

    Neff = arrayfun(@(r) r.Lin.metrics.N_eff, Rs);
    [~, imax] = max(Neff);

    Selected.(char(bs)) = Rs(imax);
end

%% =========================================================
% Plot: one row per sector, two panels
% =========================================================

secNames = fieldnames(Selected);
nsec = numel(secNames);

figure('Color','w','Position',[100 100 900 260*nsec]);

for i = 1:nsec

    secChar  = secNames{i};
    secStr   = string(secChar);
    secLabel = replace(secStr,"Sector"," sector");

    R0 = Selected.(secChar);

    % --- Get correct sector structure
    [Sstruct, Sfield] = getSectorStruct( ...
        string(R0.Sector), S_CAN, S_PAC, S_SIB, S_EUR);

    ip = R0.ProfileID;

    % --- Extract profile data
    z   = Sstruct.(Sfield).DEPTH.Value(:,ip);
    no3 = Sstruct.(Sfield).NITRAT.Value(:,ip);
    p   = Sstruct.(Sfield).PRESSURE.Value(:,ip);
    sal = Sstruct.(Sfield).SALNTY.Value(:,ip);
    wm  = Sstruct.(Sfield).watermass.Value(:,ip);
    lat = Sstruct.(Sfield).latitude.Value(ip);
    lon = Sstruct.(Sfield).longitude.Value(ip);

    % --- Convert to Absolute Salinity
    sal = gsw_SA_from_SP(sal,p,lon,lat);

    % --- Date from UNIX seconds (per profile)
    unixvec = Sstruct.(Sfield).UNIXSECS.Value(:);
    unixsec = unixvec(ip);
    
    dt = datetime(unixsec, ...
        'ConvertFrom','posixtime', ...
        'TimeZone','UTC');
    
    dateStr = datestr(dt,'dd:mmm:yyyy');

    % --- AW mask
    isAW = wm >= 30 & wm <= 50 & z >= 100 & ...
           isfinite(no3) & isfinite(sal);

    %% ---- LEFT PANEL: Nitrate vs depth
    subplot(nsec,2,2*i-1)
    plot(no3, z, 'k', 'LineWidth', 1.1); hold on
    plot(no3(isAW), z(isAW), 'r', 'LineWidth', 2)
    set(gca,'YDir','reverse')
    xlabel('Nitrate (\mumol kg^{-1})')
    ylabel('Depth (m)')
    title(secLabel + " (" + ...
          num2str(lat,'%.1f') + "°N, " + ...
          num2str(lon,'%.1f') + "°E, " + ...
          dateStr + ")")
    grid on; box on

    %% ---- RIGHT PANEL: Nitrate vs salinity (AW only)
    subplot(nsec,2,2*i)
    scatter(sal(isAW), no3(isAW), 30, 'filled'); hold on

    % --- Deming regression from Results (NO refit)
    m = R0.Deming_slope;
    b = R0.Deming_intercept;

    if isfinite(m) && isfinite(b) && nnz(isAW) >= 2
        xs = linspace(min(sal(isAW)), max(sal(isAW)), 50);
        ys = m*xs + b;
        plot(xs, ys, 'r', 'LineWidth', 2)
        slopeStr = "Deming slope = " + num2str(m,'%.2f');
    else
        slopeStr = "Deming slope = n/a";
    end

    xlabel('SA')
    ylabel('Nitrate (\mumol kg^{-1})')
    grid on; box on
    text(0.05,0.9, slopeStr, 'Units','normalized')

end

%% =========================================================
% Helper function (STRING-SAFE)
% =========================================================
function [Sstruct, Sfield] = getSectorStruct( ...
    sectorName, S_CAN, S_PAC, S_SIB, S_EUR)

    sectorName = string(sectorName);

    if contains(sectorName,"European")
        Sstruct = S_EUR;
        Sfield  = 'EuropeanSector';

    elseif contains(sectorName,"Barents")
        Sstruct = S_EUR;
        Sfield  = 'BarentsSector';

    elseif contains(sectorName,"Canadian")
        Sstruct = S_CAN;
        Sfield  = 'CanadianSector8090N';

    elseif contains(sectorName,"Pacific")
        Sstruct = S_PAC;
        Sfield  = 'PacificSector8090N';

    elseif contains(sectorName,"Siberian")
        Sstruct = S_SIB;
        Sfield  = 'SiberianSector8090N';

    else
        error("Unknown sector: " + sectorName)
    end
end