function s = S_eval(coeffs, N, pt)
% Fourier series evaluation function
% Inputs:
%  coeffs    - [ (2*N1+1) x ... x (2*Nn+1) ] array of Fourier coefficients
%  N         - Row vector [1 x n] of maximum frequency indices (|k| <= N(i))
%  pt        - points needed
% Outputs:
%  s         - the function value at pt
    s = coeffs;
    n = length(pt);
    for dim = 1:n
        kj = -N(dim):N(dim);
        phase = 1i * kj * pt(dim);
        v = exp(phase(:));
        vsize = ones(1, n);
        vsize(dim) = length(v);
        v = reshape(v, vsize);
        s = sum(s .* v, dim);
    end
    s = real(s); % Optional: use if f is real-valued
end

