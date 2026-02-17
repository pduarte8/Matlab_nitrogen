function out = diagnose_linearity_adaptive(x,y)
% DIAGNOSE_LINEARITY_ADAPTIVE
% Adaptive linear diagnostics
%
% Rules:
% - Correlation significance ALWAYS evaluated using N_eff = min(N,50)
% - Small N (<20): correlation + residual–predictor correlation
% - Medium N (20–49): correlation only
% - Large N (>=50): correlation + bootstrap slope stability
%
% Curvature is NOT handled here

N = numel(x);

% Initialize output
out = struct();
out.is_linear = false;
out.confidence = 'low';
out.reason = '';
out.metrics = struct();

% Pearson correlation
r = corr(x,y,'rows','complete');

% -------------------------
% Effective sample size for correlation test
% -------------------------
N_eff = min(N,50);

df = N_eff - 2;
tcrit = tinv(0.975, df);                 % two-sided alpha = 0.05
rcrit = tcrit / sqrt(tcrit^2 + df);

out.metrics.r = r;
out.metrics.rcrit = rcrit;
out.metrics.N_eff = N_eff;

corr_significant = abs(r) >= rcrit;

% -------------------------
% Regime definition
% -------------------------
if N >= 50
    regime = 'large';
elseif N >= 20
    regime = 'medium';
else
    regime = 'small';
end

% -------------------------
% Bootstrap parameters (large N only)
% -------------------------
boot_frac = 0.2;
boot_min  = 20;
n_boot    = 200;
lambda    = 1;

switch regime

    case 'large'
        % ======================
        % Large sample
        % ======================
        n_sub = max(round(boot_frac * N), boot_min);
        slopes = NaN(n_boot,1);

        for k = 1:n_boot
            ii = randperm(N,n_sub);
            xs = x(ii);
            ys = y(ii);

            S = cov(xs,ys,1);
            Sxx = S(1,1);
            Syy = S(2,2);
            Sxy = S(1,2);

            slopes(k) = (Syy - lambda*Sxx + ...
                sqrt((Syy - lambda*Sxx)^2 + 4*lambda*Sxy^2)) ...
                / (2*Sxy);
        end

        slopes = slopes(isfinite(slopes));
        beta_mean = mean(slopes);
        beta_cv   = std(slopes) / abs(beta_mean);
        
        sign_consistency = mean(sign(slopes) == sign(beta_mean));
        
        out.metrics.beta_cv = beta_cv;
        out.metrics.sign_consistency = sign_consistency;
        
        if corr_significant && beta_cv < 0.25 && sign_consistency >= 0.9
            out.is_linear = true;
            out.confidence = 'high';
            out.reason = 'large-N: significant correlation, stable and sign-consistent bootstrap slopes';
        else
            out.reason = 'large-N: unstable or sign-inconsistent bootstrap slopes';
        end

    case 'medium'
        % ======================
        % Medium sample
        % ======================
        if corr_significant
            % Linear fit
            p = polyfit(x,y,1);
            yhat = polyval(p,x);
            res = y - yhat;

            % Residual–predictor correlation
            r_res = corr(res,x,'rows','complete');
            out.metrics.r_res = r_res;

            if abs(r_res) < rcrit
                out.is_linear = true;
                out.confidence = 'low';
                out.reason = 'medium-N: significant correlation and uncorrelated residuals';
            else
                out.reason = 'medium-N: residuals correlated with predictor';
            end
        else
            out.reason = 'medium-N: correlation not significant';
        end
     
    case 'small'
        % ======================
        % Small sample
        % ======================
        if corr_significant
            % Linear fit
            p = polyfit(x,y,1);
            yhat = polyval(p,x);
            res = y - yhat;

            % Residual–predictor correlation
            r_res = corr(res,x,'rows','complete');
            out.metrics.r_res = r_res;

            if abs(r_res) < rcrit
                out.is_linear = true;
                out.confidence = 'low';
                out.reason = 'small-N: significant correlation and uncorrelated residuals';
            else
                out.reason = 'small-N: residuals correlated with predictor';
            end
        else
            out.reason = 'small-N: correlation not significant';
        end
end

end