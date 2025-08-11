function coeffs = MultiFourierSeries(f, a, b, N, M, quad_method)
% MULTIFOURIERSERIES Truncated Fourier series for multivariable functions using quadratures.
%   COEFFS = MultiFourierSeries(F, A, B, N, M, QUAD_METHOD, ORDER) computes Fourier coefficients
%   for function F on the hypercube [A,B]^n. N sets frequency limits, M sets integration points,
%   QUAD_METHOD chooses quadrature type ('newton_cotes', 'gauss_legendre', 'clenshaw_curtis'),
%   and ORDER specifies Newton-Cotes order (1 for trapezoidal, 2 for Simpson).
%
% Inputs:
%   f          - Function handle (accepts [P x n] matrix, returns [P x 1] vector)
%   a, b       - Row vectors [1 x n] of lower/upper bounds
%   N          - Row vector [1 x n] of maximum frequency indices
%   M          - (Optional) Row vector [1 x n] of integration points per dimension
%   quad_method- (Optional) Quadrature method: 'newton_cotes' (default), 'gauss_legendre', 'clenshaw_curtis'
%   order      - (Optional) Newton-Cotes order (default=2)
%
% Outputs:
%   coeffs - [ (2*N1+1) x ... x (2*Nn+1) ] array of Fourier coefficients

n = length(a); % Number of dimensions
V = prod(b - a); % Hypercube volume

% Default integration points
if nargin < 5 || isempty(M)
    M = 10 * (2*N + 1);
end

% Default quadrature method and order
if nargin < 6 || isempty(quad_method)
    quad_method = 'nc';
end
% Generate quadrature points and weights for each dimension
x = cell(1, n);
w1D = cell(1, n);
for j = 1:n
    if strcmp(quad_method, 'gl')
        [x{j}, w1D{j}] = legpts(M(j), [a(j), b(j)]);
        w1D{j} = w1D{j}';
    elseif strcmp(quad_method, 'cc')
        [x{j}, w1D{j}] = chebpts(M(j), [a(j), b(j)]);
        w1D{j} = w1D{j}';
    elseif strcmp(quad_method, 'nc')
        [x{j}, w1D{j}] = trigpts(M(j), [a(j), b(j)]);
        w1D{j} = w1D{j}';
    elseif strcmp(quad_method, 'gjl')
        [x{j}, w1D{j}] = lobpts(M(j));
        x{j} = (b(j) - a(j)) * (x{j} + 1)/2 + a(j);
        w1D{j} = w1D{j}';
    end
end

% Build grid and evaluate function
[Xgrid{1:n}] = ndgrid(x{:});
points = reshape(cat(n+1, Xgrid{:}), [], n);
Fvals = f(points);
Fvals = reshape(Fvals, M); % Reshape to grid dimensions

% Construct weight tensor
W = ones(M);
for j = 1:n
    w_shape = ones(1, n);
    w_shape(j) = M(j);
    Wj = reshape(w1D{j}, w_shape);
    W = W .* Wj; % Tensor product of weights
end
Fvals_weighted = Fvals .* W; % Apply weights

% Compute Fourier coefficients via dimension-wise transforms
coeffs = Fvals_weighted;
for j = 1:n
    kj = (-N(j):N(j))';
    E = exp(-1i * kj * x{j}'); % Basis functions
    
    % Move dimension j to front for transformation
    perm = [j, 1:j-1, j+1:n];
    coeffs = permute(coeffs, perm);
    siz = size(coeffs);
    coeffs = reshape(coeffs, siz(1), []);
    coeffs = E * coeffs; % Transform along dimension j
    siz(1) = length(kj);
    coeffs = reshape(coeffs, siz);
end

% Reorder dimensions to match original variable order
if n > 1
    coeffs = permute(coeffs, n:-1:1);
end

% Normalize by hypercube volume
coeffs = coeffs / V;
end