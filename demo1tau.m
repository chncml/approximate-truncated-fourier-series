function demo1tau
%% parameter choices in our algorithms
   clc;clear;
   epsilon = 1e-7;  % Tolerance
   kmax = 100;       % Max iterations

   alpha = 2;  % Smoothness parameter (assumed known)
   C_constant = 1; % Constant from Theorem 3.1
   I1 = ceil((C_constant / epsilon)^(1/alpha));
   I2 = I1;

   a = [-pi, -pi];b = [pi, pi];

   % Truncation: 5 terms per dimension (|k_x|<=2, |k_y|<=2)
   N = [I1, I2];

   % Integration resolution
   M = [5000, 5000] + 1; % 100 points per dimension

   quad_method = 'cc';
   b1 = 4;b2 = 4;

   % quad_method = 'gl';
   % b1 = 6;b2 = 6;
   % 
   % quad_method = 'nc';
   % b1 = 6;b2 = 6;
   
   % the points
   K1 = 60;x = linspace(-pi, pi, K1);
   K2 = 60;y = linspace(-pi, pi, K2);
   [X, Y] = ndgrid(x, y);

   %% target function
   g = @(p) (5^(3/4) * 15/4/sqrt(3))^2 * max(0, 1/5 - (p(:, 1) - 1/2).^2) .* max(0, 1/5 - (p(:, 2) - 1/2).^2);
   f = @(x, y) (5^(3/4) * 15/4/sqrt(3))^2 * max(0, 1/5 - (x - 1/2).^2) .* max(0, 1/5 - (y - 1/2).^2);

   f_vals = f(X, Y);
   %% tau=1e-1
   t1 = tic;
   [C1, U1, R1] = adaptiveFourierApproximation_quad(f, b1, b2, 1e-1, kmax, I1, I2, M(1), M(2), quad_method);
   coeffs1 = C1 * (U1 * R1);
   t1 = toc(t1);

   f_approx1 = zeros(K1, K2);
   for i1 = 1 : K1
       i1
       for i2 = 1 : K2
           i2
           f_approx1(i1, i2) = S_eval(coeffs1, N, [X(i1, i2), Y(i1, i2)]);
       end
   end

   %% tau=1e-2
   t2 = tic;
   [C2, U2, R2] = adaptiveFourierApproximation_quad(f, b1, b2, 1e-2, kmax, I1, I2, M(1), M(2), quad_method);
   coeffs2 = C2 * (U2 * R2);
   t2 = toc(t2);

   f_approx2 = zeros(K1, K2);
   for i1 = 1 : K1
       i1
       for i2 = 1 : K2
           i2
           f_approx2(i1, i2) = S_eval(coeffs2, N, [X(i1, i2), Y(i1, i2)]);
       end
   end

   %% tau=1e-3
   t3 = tic;
   [C3, U3, R3] = adaptiveFourierApproximation_quad(f, b1, b2, 1e-3, kmax, I1, I2, M(1), M(2), quad_method);
   coeffs3 = C3 * (U3 * R3);
   t3 = toc(t3);

   f_approx3 = zeros(K1, K2);
   for i1 = 1 : K1
       i1
       for i2 = 1 : K2
           i2
           f_approx3(i1, i2) = S_eval(coeffs3, N, [X(i1, i2), Y(i1, i2)]);
       end
   end

   %%  tau=1e-4
   t4 = tic;
   [C4, U4, R4] = adaptiveFourierApproximation_quad(f, b1, b2, 1e-4, kmax, I1, I2, M(1), M(2), quad_method);
   coeffs4 = C4 * (U4 * R4);
   t4 = toc(t4);

   f_approx4 = zeros(K1, K2);
   for i1 = 1 : K1
       i1
       for i2 = 1 : K2
           i2
           f_approx4(i1, i2) = S_eval(coeffs4, N, [X(i1, i2), Y(i1, i2)]);
       end
   end

   %% tau=1e-5
   t5 = tic;
   [C5, U5, R5] = adaptiveFourierApproximation_quad(f, b1, b2, 1e-5, kmax, I1, I2, M(1), M(2), quad_method);
   coeffs5 = C5 * (U5 * R5);
   t5 = toc(t5);

   f_approx5 = zeros(K1, K2);
   for i1 = 1 : K1
       i1
       for i2 = 1 : K2
           i2
           f_approx5(i1, i2) = S_eval(coeffs5, N, [X(i1, i2), Y(i1, i2)]);
       end
   end

   %% tau=1e-5
   t6 = tic;
   [C6, U6, R6] = adaptiveFourierApproximation_quad(f, b1, b2, 1e-6, kmax, I1, I2, M(1), M(2), quad_method);
   coeffs6 = C6 * (U6 * R6);
   t6 = toc(t6);

   f_approx6 = zeros(K1, K2);
   for i1 = 1 : K1
       i1
       for i2 = 1 : K2
           i2
           f_approx6(i1, i2) = S_eval(coeffs6, N, [X(i1, i2), Y(i1, i2)]);
       end
   end

   %% tau=1e-7
   t7 = tic;
   [C7, U7, R7] = adaptiveFourierApproximation_quad(f, b1, b2, 1e-7, kmax, I1, I2, M(1), M(2), quad_method);
   coeffs7 = C7 * (U7 * R7);
   t7 = toc(t7);

   f_approx7 = zeros(K1, K2);
   for i1 = 1 : K1
       i1
       for i2 = 1 : K2
           i2
           f_approx7(i1, i2) = S_eval(coeffs7, N, [X(i1, i2), Y(i1, i2)]);
       end
   end

   %% tau=1e-8
   t8 = tic;
   [C8, U8, R8] = adaptiveFourierApproximation_quad(f, b1, b2, 1e-8, kmax, I1, I2, M(1), M(2), quad_method);
   coeffs8 = C8 * (U8 * R8);
   t8 = toc(t8);

   f_approx8 = zeros(K1, K2);
   for i1 = 1 : K1
       i1
       for i2 = 1 : K2
           i2
           f_approx8(i1, i2) = S_eval(coeffs8, N, [X(i1, i2), Y(i1, i2)]);
       end
   end

   %% tau=1e-9
   t9 = tic;
   [C9, U9, R9] = adaptiveFourierApproximation_quad(f, b1, b2, 1e-9, kmax, I1, I2, M(1), M(2), quad_method);
   coeffs9 = C9 * (U9 * R9);
   t9 = toc(t9);

   f_approx9 = zeros(K1, K2);
   for i1 = 1 : K1
       i1
       for i2 = 1 : K2
           i2
           f_approx9(i1, i2) = S_eval(coeffs9, N, [X(i1, i2), Y(i1, i2)]);
       end
   end

   %% tau=1e-10
   t10 = tic;
   [C10, U10, R10] = adaptiveFourierApproximation_quad(f, b1, b2, 1e-10, kmax, I1, I2, M(1), M(2), quad_method);
   coeffs10 = C10 * (U10 * R10);
   t10 = toc(t10);

   f_approx10 = zeros(K1, K2);
   for i1 = 1 : K1
       i1
       for i2 = 1 : K2
           i2
           f_approx10(i1, i2) = S_eval(coeffs10, N, [X(i1, i2), Y(i1, i2)]);
       end
   end

   subplot(2, 5, 1)
   mesh(X, Y, abs(f_approx1 - f_vals));
   hold on
   subplot(2, 5, 2)
   mesh(X, Y, abs(f_approx2 - f_vals));
   hold on
   subplot(2, 5, 3)
   mesh(X, Y, abs(f_approx3 - f_vals));
   hold on
   subplot(2, 5, 4)
   mesh(X, Y, abs(f_approx4 - f_vals));
   hold on
   subplot(2, 5, 5)
   mesh(X, Y, abs(f_approx5 - f_vals));
   hold on
   subplot(2, 5, 6)
   mesh(X, Y, abs(f_approx6 - f_vals));
   hold on
   subplot(2, 5, 7)
   mesh(X, Y, abs(f_approx7 - f_vals));
   hold on
   subplot(2, 5, 8)
   mesh(X, Y, abs(f_approx8 - f_vals));
   hold on
   subplot(2, 5, 9)
   mesh(X, Y, abs(f_approx9 - f_vals));
   hold on
   subplot(2, 5, 10)
   mesh(X, Y, abs(f_approx10 - f_vals));
   hold on

   [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10]
end