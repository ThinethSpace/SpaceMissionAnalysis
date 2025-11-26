classdef functionsOrbit
    methods(Static)
        %% 1a) Gauss-Lambert (Universal-Variable-Form, numerische Ableitung)
        function [v1, v2] = lambertGauss(r1, r2, dt, mu)
            r1n = norm(r1);
            r2n = norm(r2);
            cos_dtheta = dot(r1, r2) / (r1n * r2n);
            dtheta = acos(cos_dtheta);

            % Kurzweg-Transfer
            A = sin(dtheta) * sqrt(r1n * r2n / (1 - cos_dtheta));

            % Startwert für z
            z = 0;
            tol = 1e-8;
            dz_num = 1e-5;

            for k = 1:100
                [C, S] = functionsOrbit.stumpff(z);
                y = r1n + r2n + A * ((z * S - 1) / sqrt(C));
                F = ( (y / C)^(1.5) ) * S + A * sqrt(y) - sqrt(mu) * dt;

                if abs(F) < tol
                    break
                end

                % numerische Ableitung dF/dz
                z_plus = z + dz_num;
                z_minus = z - dz_num;

                [C_p, S_p] = functionsOrbit.stumpff(z_plus);
                y_p = r1n + r2n + A * ((z_plus * S_p - 1) / sqrt(C_p));
                F_p = ( (y_p / C_p)^(1.5) ) * S_p + A * sqrt(y_p) - sqrt(mu) * dt;

                [C_m, S_m] = functionsOrbit.stumpff(z_minus);
                y_m = r1n + r2n + A * ((z_minus * S_m - 1) / sqrt(C_m));
                F_m = ( (y_m / C_m)^(1.5) ) * S_m + A * sqrt(y_m) - sqrt(mu) * dt;

                dFdz = (F_p - F_m) / (2 * dz_num);

                z = z - F / dFdz;
            end

            % Endgültige Größen
            [C, S] = functionsOrbit.stumpff(z);
            y = r1n + r2n + A * ((z * S - 1) / sqrt(C));

            f = 1 - y / r1n;
            g = A * sqrt(y / mu);
            gdot = 1 - y / r2n;

            v1 = (r2 - f * r1) / g;
            v2 = (gdot * r2 - r1) / g;  % v2 ist hier die Ankunftsgeschwindigkeit (in der Aufgabe v3 genannt)
        end

        %% Stumpff-Funktionen C(z), S(z)
        function [C, S] = stumpff(z)
            if z > 0
                s = sqrt(z);
                C = (1 - cos(s)) / z;
                S = (s - sin(s)) / (s^3);
            elseif z < 0
                s = sqrt(-z);
                C = (1 - cosh(s)) / z;  % z < 0
                S = (sinh(s) - s) / (s^3);
            else
                C = 0.5;
                S = 1/6;
            end
        end

        %% RV -> Kepler-Elemente
        function [a, e, i, RAAN, omega, nu] = rv2coe(r, v, mu)
            rnorm = norm(r);
            vnorm = norm(v);

            h = cross(r, v);
            hnorm = norm(h);

            kvec = [0; 0; 1];
            i = acos(h(3) / hnorm);

            n = cross(kvec, h);
            nnorm = norm(n);

            evec = (1/mu) * ( (vnorm^2 - mu/rnorm) * r - dot(r, v) * v );
            e = norm(evec);

            RAAN = acos(n(1) / nnorm);
            if n(2) < 0
                RAAN = 2*pi - RAAN;
            end

            omega = acos(dot(n, evec) / (nnorm * e));
            if evec(3) < 0
                omega = 2*pi - omega;
            end

            nu = acos(dot(evec, r) / (e * rnorm));
            if dot(r, v) < 0
                nu = 2*pi - nu;
            end

            a = 1 / (2/rnorm - vnorm^2/mu);
        end

        %% 4th-Order Runge-Kutta für 2-Körper
        function [tArr, rArr, vArr] = propagateRK4(r0, v0, mu, dtTotal, dtStep)
            nSteps = floor(dtTotal / dtStep) + 1;
            tArr  = linspace(0, dtTotal, nSteps);
            rArr  = zeros(3, nSteps);
            vArr  = zeros(3, nSteps);

            rArr(:,1) = r0;
            vArr(:,1) = v0;

            for k = 1:(nSteps-1)
                y = [rArr(:,k); vArr(:,k)];

                k1 = functionsOrbit.twoBodyRHS(y, mu);
                k2 = functionsOrbit.twoBodyRHS(y + 0.5*dtStep*k1, mu);
                k3 = functionsOrbit.twoBodyRHS(y + 0.5*dtStep*k2, mu);
                k4 = functionsOrbit.twoBodyRHS(y + dtStep*k3, mu);

                y_next = y + dtStep/6 * (k1 + 2*k2 + 2*k3 + k4);

                rArr(:,k+1) = y_next(1:3);
                vArr(:,k+1) = y_next(4:6);
            end
        end

        %% RHS des 2-Körper-Problems
        function dydt = twoBodyRHS(y, mu)
            r = y(1:3);
            v = y(4:6);
            rnorm = norm(r);

            a = -mu * r / (rnorm^3);

            dydt = [v; a];
        end

        %% ==================== TASK 2: Gibbs-Methode =======================
        function v2 = gibbs(r1, r2, r3, mu)
            % Gibbs-Methode zur Bestimmung von v2
            % r1, r2, r3 : Positionsvektoren (3x1) [km]
            % mu         : Gravitationsparameter [km^3/s^2]
            %
            % Rückgabe:
            % v2         : Geschwindigkeitsvektor im zweiten Punkt [km/s]

            % Kreuzprodukte
            c12 = cross(r1, r2);
            c23 = cross(r2, r3);
            c31 = cross(r3, r1);

            % D-, N- und S-Vektoren nach klassischer Gibbs-Formel
            D = c12 + c23 + c31;

            N = norm(r1)*c23 + norm(r2)*c31 + norm(r3)*c12;

            S = r1*(norm(r2) - norm(r3)) + ...
                r2*(norm(r3) - norm(r1)) + ...
                r3*(norm(r1) - norm(r2));

            % v2 aus Gibbs
            v2 = sqrt(mu/(norm(N)*norm(D))) * ( cross(D, r2)/norm(r2) + S );
        end

        %% ==================== TASK 3: Herrick-Gibbs =======================
        function v2 = herrickGibbs(r1, r2, r3, t1, t2, t3, mu)
            % Herrick-Gibbs-Methode zur Bestimmung von v2 (am mittleren Zeitpunkt)
            % r1, r2, r3 : Positionsvektoren (3x1) [km]
            % t1, t2, t3 : Zeitpunkte [s]
            % mu         : Gravitationsparameter [km^3/s^2]
            %
            % Rückgabe:
            % v2         : Geschwindigkeitsvektor im zweiten Punkt [km/s]

            % Zeitabstände
            tau21 = t2 - t1;    % t2 - t1
            tau31 = t3 - t1;    % t3 - t1
            tau32 = t3 - t2;    % t3 - t2

            muOTwelve = mu / 12.0;

            r1n = norm(r1);
            r2n = norm(r2);
            r3n = norm(r3);

            % Koeffizienten nach Herrick-Gibbs (vgl. Vallado / Orekit) :contentReference[oaicite:0]{index=0}
            c1 = -tau32 * (1/(tau21 * tau31) + muOTwelve / (r1n^3));
            c2 = (tau32 - tau21) * (1/(tau21 * tau32) + muOTwelve / (r2n^3));
            c3 =  tau21 * (1/(tau32 * tau31) + muOTwelve / (r3n^3));

            % v2 als Linearkombination der Positionsvektoren
            v2 = c1 * r1 + c2 * r2 + c3 * r3;
        end
    end
end
