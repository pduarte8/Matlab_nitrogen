
function [slope, intercept] = deming_regress(x, y, lambda)
    if nargin < 3, lambda = 1; end  % ratio of variances (σ_y²/σ_x²), =1 when unknown
    x = x(:); y = y(:);
    n = numel(x);
    mx = mean(x); my = mean(y);
    Sxx = sum((x-mx).^2)/(n-1);
    Syy = sum((y-my).^2)/(n-1);
    Sxy = sum((x-mx).*(y-my))/(n-1);

    slope = (Syy - lambda*Sxx + sqrt((Syy - lambda*Sxx).^2 + 4*lambda*Sxy.^2)) / (2*Sxy);
    intercept = my - slope*mx;
end