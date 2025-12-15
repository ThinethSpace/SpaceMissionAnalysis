% Test secant method
% Supergolden ratio test

IT = InterplanetaryTransfers;

f = @(x) x^3 - x^2 -1;

[x_sol, x_all] = IT.perform_secant_method(f, 1, 2, 20, 1E-10);

