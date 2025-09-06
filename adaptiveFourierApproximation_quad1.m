function [C, U, R] = adaptiveFourierApproximation_quad1(f, b1, b2, tau, kmax, I1, I2, M1, M2, quad_method)
    % Inputs:
    %   f: Function handle for f(x1, x2)
    %   b1, b2: Block sizes for mode 1 and 2
    %   tau: Tolerance for stopping criterion
    %   kmax: Maximum number of iterations
    %   I1, I2: Maximum frequency indices for each mode
    %   M1, M2: (Optional) Number of quadrature points for each mode
    %   quad_type: Quadrature type ('NC', 'GL', 'CC')
    %
    % Outputs:
    %   C, U, R: Approximation matrices
    %
    % Handle optional inputs
    if nargin < 8
        M1 = 2*I1 + 1;
        M2 = 2*I2 + 1;
    end
    if nargin < 10 || isempty(quad_method)
        quad_method = 'nc'; % Default to Newton-Cotes
    end

    % Generate quadrature points and weights
    if strcmp(quad_method, 'gl')
        [x1, w1] = legpts(M1, [-pi, pi]);w1 = w1';
        [x2, w2] = legpts(M2, [-pi, pi]);w2 = w2';
    elseif strcmp(quad_method, 'cc')
        [x1, w1] = chebpts(M1, [-pi, pi]);w1 = w1';
        [x2, w2] = chebpts(M2, [-pi, pi]);w2 = w2';
    elseif strcmp(quad_method, 'nc')
        [x1, w1] = trigpts(M1, [-pi, pi]);w1 = w1';
        [x2, w2] = trigpts(M2, [-pi, pi]);w2 = w2';
    elseif strcmp(quad_method, 'gjl')
        [x1, w1] = lobpts(M1);x1 = pi * (x1 + 1) - pi;w1 = w1';
        [x2, w2] = lobpts(M2);x2 = pi * (x2 + 1) - pi;w2 = w2';
    end
    
    % Evaluate function on the grid
    [X1, X2] = ndgrid(x1, x2);
    F_values = f(X1, X2);
    
    % Apply quadrature weights
    F_weighted = F_values .* (w1 * w2');
    
    % All frequency indices
    k1_all = (-I1:I1)';
    k2_all = (-I2:I2)';
    
    % Helper function to compute Fourier coefficients
    getCoeffs = @(k1_vals, k2_vals) (1/(4*pi^2)) * (exp(-1i * k1_vals * x1') * F_weighted * exp(-1i * x2 * k2_vals'));
    
    % Initialize index sets and core matrix
    I1_current = [];
    I2_current = [];
    C = [];R = [];
    nF_C = 0;nF_R = 0;
    b1half = ceil(b1/2);b2half = ceil(b2/2);
    
    % Randomly select initial indices
    % idx1 = randperm(length(k1_all), min(b1, length(k1_all)));
    % idx2 = randperm(length(k2_all), min(b2, length(k2_all)));
    idx1 = [((I1 + 1) - b1half) : ((I1 + 1) + b1half)];
    idx2 = [((I2 + 1) - b2half) : ((I2 + 1) + b2half)];
    I1_current = k1_all(idx1);
    I2_current = k2_all(idx2);
    
    % Compute initial C and R
    C = getCoeffs(k1_all, I2_current);
    R = getCoeffs(I1_current, k2_all);
    nF_C = norm(C, 'fro')^2;nF_R = norm(R, 'fro')^2;
    sc = svd(C, 'econ');sr = svd(R, 'econ');
    tol = max(min(sc) / sqrt(nF_C), min(sr) / sqrt(nF_R));
    
    % Iteratively add indices
    for iter = 1:kmax
        if tol <= tau
            break;
        end
        
        % Find available frequencies
        avail_k1 = setdiff(k1_all, I1_current);
        avail_k2 = setdiff(k2_all, I2_current);
        if isempty(avail_k1) || isempty(avail_k2)
            break;
        end
        
        % Select new indices
        % idx1 = randperm(length(avail_k1), min(b1, length(avail_k1)));
        % idx2 = randperm(length(avail_k2), min(b2, length(avail_k2)));
        idx1 = [(I1 + 1 - (iter + 1) * b1half) : (I1 + 1 - iter * b1half - 1), (I1 + 1 + iter * b1half + 1) : (I1 + 1 + (iter + 1) * b1half)];
        idx2 = [(I2 + 1 - (iter + 1) * b2half) : (I2 + 1 - iter * b2half - 1), (I2 + 1 + iter * b2half + 1) : (I2 + 1 + (iter + 1) * b2half)];
        new_k1 = avail_k1(idx1);
        new_k2 = avail_k2(idx2);
        
        % Compute new blocks
        C1k = getCoeffs(k1_all, new_k2);
        R1k = getCoeffs(new_k2, k2_all);
        
        % Update core matrix and Frobenius norm
        nF_C = nF_C + norm(C1k, 'fro')^2;
        nF_R = nF_R + norm(R1k, 'fro')^2;
        C = [C, C1k];R = [R; R1k];
        
        % Update index sets
        I1_current = [I1_current; new_k1];
        I2_current = [I2_current; new_k2];
        
        % Update tolerance
        sc = svd(C, 'econ');sr = svd(R, 'econ');
        tol = max(min(sc) / sqrt(nF_C), min(sr) / sqrt(nF_R));% Update tolerance
    end

    % Compute pseudoinverse of core matrix
    % G = getCoeffs(I1_current, I2_current);
    G = C(I1_current + I1 +1, :);
    U = pinv(G);
end
