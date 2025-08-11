function demo2comparison
%% parameter setting
   clc;clear;
   epsilon = 1e-7;  % Tolerance
   tau = 1e-5;      % Stopping tolerance
   kmax = 100;       % Max iterations

   alpha = 2;  % Smoothness parameter (assumed known)
   C_constant = 1; % Constant from Theorem 3.1
   I1 = ceil((C_constant / epsilon)^(1/alpha));
   I2 = I1;
   [I1, I2]

   a = [-pi, -pi];b = [pi, pi];

   % Truncation: 5 terms per dimension (|k_x|<=2, |k_y|<=2)
   N = [I1, I2];

   % Integration resolution
   M = [5000, 5000] + 1; % 100 points per dimension

   quad_method = 'cc';
   b1 = 6;b2 = 6;

   % quad_method = 'gl';
   % b1 = 14;b2 = 14;
   % 
   % quad_method = 'nc';
   % b1 = 14;b2 = 14;
   % the points
   K1 = 60;x = linspace(-pi, pi, K1);
   K2 = 60;y = linspace(-pi, pi, K2);
   [X, Y] = ndgrid(x, y);

   %% target function
   g = @(p) ((1 - p(:,1).^2 - p(:,2).^2) .* exp(p(:,1).* cos(p(:,2))));
   f = @(x1, x2) ((1 - x1.^2 - x2.^2) .* exp(x1.* cos(x2)));
   f_vals = f(X, Y);

   %% truncated fourier series
   T1 = zeros(10, 1);Coeffs1 = zeros(2 * I1 + 1, 2 * I2 + 1);
   for sample = 1 : 10
       t = tic;coeffs1 = MultiFourierSeries(g, a, b, N, M, quad_method);T1(sample) = toc(t);
       Coeffs1 = Coeffs1 + coeffs1;
   end
   Coeffs1 = Coeffs1/10;
   f_approx1 = zeros(K1, K2);
   for i1 = 1 : K1
       i1
       for i2 = 1 : K2
           i2
           f_approx1(i1, i2) = S_eval(Coeffs1, N, [X(i1, i2), Y(i1, i2)]);
       end
   end

   %% the first approximation method
   T2 = zeros(10, 1);Coeffs2 = zeros(2 * I1 + 1, 2 * I2 + 1);
   for sample = 1 : 10
       t = tic;
       [C2, U2, R2] = adaptiveFourierApproximation_quad1(f, b1, b2, tau, kmax, I1, I2, M(1), M(2), quad_method);
       coeffs2 = C2 * (U2 * R2);
       T2(sample) = toc(t);
       Coeffs2 = Coeffs2 + coeffs2;
   end
   Coeffs2 = Coeffs2/10;
   f_approx2 = zeros(K1, K2);
   for i1 = 1 : K1
       i1
       for i2 = 1 : K2
           i2
           f_approx2(i1, i2) = S_eval(Coeffs2, N, [X(i1, i2), Y(i1, i2)]);
       end
   end

   %% the second approximation method
   T3 = zeros(10, 1);Coeffs3 = zeros(2 * I1 + 1, 2 * I2 + 1);
   for sample = 1 : 10
       t = tic;
       [C3, U3, R3] = adaptiveFourierApproximation_quad(f, b1, b2, tau, kmax, I1, I2, M(1), M(2), quad_method);
       coeffs3 = C3 * (U3 * R3);
       T3(sample) = toc(t);
       Coeffs3 = Coeffs3 + coeffs3;
   end
   Coeffs3 = Coeffs3/10;
   f_approx3 = zeros(K1, K2);
   for i1 = 1 : K1
       i1
       for i2 = 1 : K2
           i2
           f_approx3(i1, i2) = S_eval(Coeffs3, N, [X(i1, i2), Y(i1, i2)]);
       end
   end

   subplot(241)
   mesh(X, Y, f_vals);
   hold on
   subplot(242)
   mesh(X, Y, f_approx1);
   hold on
   subplot(243)
   mesh(X, Y, f_approx2);
   hold on
   subplot(244)
   mesh(X, Y, f_approx3);
   hold on
   subplot(246)
   mesh(X, Y, abs(f_approx1 - f_vals));
   hold on
   subplot(247)
   mesh(X, Y, abs(f_approx2 - f_vals));
   hold on
   subplot(248)
   mesh(X, Y, abs(f_approx3 - f_vals));
   hold on

   [sum(T1), sum(T2), sum(T3)]/10
end