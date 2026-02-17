function Quad = diagnose_quadratic_structure(x, y)
% Diagnose quadratic structure for regression
% Outputs:
%   Quad.quad: intercept, linear, curvature (in original x scale)
%   Quad.passes_complementary_tests: boolean
%   Quad.metrics: additional diagnostics (vertex location, curvature magnitude)

N = numel(x);
%if N < 6
%    Quad = struct();
%    Quad.quad = struct('intercept',NaN,'linear',NaN,'curvature',NaN);
%    Quad.passes_complementary_tests = false;
%    Quad.metrics = struct();
%    return
%end

% Center x
x_mean = mean(x);
xc = x - x_mean;

% Quadratic fit on centered x
Xq = [ones(N,1) xc xc.^2];
beta_centered = Xq\y;
yhatQ = Xq*beta_centered;

% Convert coefficients back to original x scale
a = beta_centered(3); % curvature
b = beta_centered(2); % linear term on centered x
c = beta_centered(1); % intercept on centered x

intercept_orig = c - b*x_mean + a*x_mean^2;
linear_orig    = b - 2*a*x_mean;
curvature_orig = a;

Quad.quad.intercept = intercept_orig;
Quad.quad.linear    = linear_orig;
Quad.quad.curvature = curvature_orig;

% ========================
% Complementary diagnostics
% ========================
curv_mag = abs(curvature_orig);
vertex_x = -linear_orig/(2*curvature_orig);
x_range = max(x) - min(x);

% Criteria for "reasonable" quadratic:
% - curvature not extremely small (if very small, may be visually linear)
% - vertex within data range
curv_tol = 0.0; % Started with 0.1, then I tried zero and called it zero_curv_tol and found no differences in the results
% then I used 1e10 and most curvilinear regressions id not show up (I guess that remained only those that passed only AIC)
% can be tuned if needed

pass_curv = curv_mag >= curv_tol;
pass_vertex = (vertex_x >= min(x)) && (vertex_x <= max(x));

Quad.passes_complementary_tests = pass_curv && pass_vertex;

% Metrics for inspection
Quad.metrics.curvature = curv_mag;
Quad.metrics.vertex_x = vertex_x;
Quad.metrics.RSS = sum((y - yhatQ).^2);

end