% Example function: f(x) = x^2 - 2 (Root at sqrt(2))
f = @(x) x^2 - 2;
% Derivative of the function: f'(x) = 2x
df = @(x) 2*x;
% Initial guess
x0 = 1;
% Tolerance for convergence
tol = 1e-6;
% Maximum number of iterations
max_iter = 100;

%create opts struct
opts.k_max = max_iter;
opts.tol = tol;
opts.return_all = true;

% Newton's method implementation
[x1, x_all] = perform_newtons_method(f, df, x0, opts);

x1
x_all