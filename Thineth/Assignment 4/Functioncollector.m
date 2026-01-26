classdef Functioncollector
    methods (Static)

        function P = task2_defaultParams()
            % Constants (Earth)
            P.mu = 3.986004418e14;         % [m^3/s^2]
            P.Re = 6378.1363e3;            % [m] (equatorial radius is fine)

            % Initial LEO (given)
            P.h0 = 200e3;                  % [m]
            P.r0 = P.Re + P.h0;            % [m]
            P.i0 = deg2rad(22.5);          % [rad]
            P.Om = deg2rad(40);            % [rad] RAAN
            P.u0 = deg2rad(17);            % [rad] argument of latitude at t0

            % Target GEO + initial longitude
            P.rGEO = 42164e3;              % [m] standard GEO radius
            P.lambda_t0 = deg2rad(200);    % [rad]

            % GEO mean motion
            P.nGEO = sqrt(P.mu / P.rGEO^3);

            % LEO mean motion (circular)
            P.nLEO = sqrt(P.mu / P.r0^3);

            % Transfer ellipse parameters (standard LEO->GEO Hohmann-like)
            P.aT   = 0.5*(P.r0 + P.rGEO);
            P.tTOF = pi*sqrt(P.aT^3 / P.mu);

            % Speeds used by assignment at burn points
            P.vLEO = sqrt(P.mu / P.r0);                               % circular at r0
            P.vGEO = sqrt(P.mu / P.rGEO);                             % circular at rGEO
            P.vPer = sqrt(P.mu*(2/P.r0   - 1/P.aT));                  % transfer at perigee
            P.vApo = sqrt(P.mu*(2/P.rGEO - 1/P.aT));                  % transfer at apogee

            % Longitude mismatch weight: choose so that ~1 deg matters ~O(1..10 m/s)
            % penalty = w * (dLambda_rad)^2
            % For 1 deg = 0.01745 rad -> squared ~3.046e-4
            % w=5000 => penalty ~15.2 (comparable to tens of m/s scale)
            P.wlambda = 5000;
        end


        function J = task2_objective(z, P)
            % z = [tw; eta]
            tw  = z(1);
            eta = z(2);

            % enforce bounds softly if called outside fmincon (robustness)
            if tw < 0, J = 1e30; return; end
            if eta < 0 || eta > 1, J = 1e30; return; end

            % total DV from spec (depends on eta only, with the provided speed model)
            dI1 = eta * P.i0;
            dI2 = (1 - eta) * P.i0;

            dV1 = Functioncollector.dv_combined(P.vLEO, P.vPer, dI1);
            dV2 = Functioncollector.dv_combined(P.vApo, P.vGEO, dI2);
            dVtot = dV1 + dV2;

            % longitude mismatch at arrival
            tArr = tw + P.tTOF;

            lambda_t = P.lambda_t0 + P.nGEO*tArr;        % target
            lambda_s = Functioncollector.interceptor_arrival_longitude(tw, eta, P);

            dLambda = Functioncollector.wrapToPi(lambda_s - lambda_t);

            % objective
            J = dVtot + P.wlambda*(dLambda^2);
        end


        function R = task2_evaluate(z, P)
            % Convenience: compute all reportable outputs
            R.tw  = z(1);
            R.eta = z(2);

            R.tTOF = P.tTOF;
            R.ttot = R.tw + R.tTOF;

            dI1 = R.eta * P.i0;
            dI2 = (1 - R.eta) * P.i0;

            R.dV1 = Functioncollector.dv_combined(P.vLEO, P.vPer, dI1);
            R.dV2 = Functioncollector.dv_combined(P.vApo, P.vGEO, dI2);
            R.dVtot = R.dV1 + R.dV2;

            lambda_t = P.lambda_t0 + P.nGEO*R.ttot;
            lambda_s = Functioncollector.interceptor_arrival_longitude(R.tw, R.eta, P);

            R.lambda_t = Functioncollector.wrapToPi(lambda_t);
            R.lambda_s = Functioncollector.wrapToPi(lambda_s);

            dLambda = Functioncollector.wrapToPi(R.lambda_s - R.lambda_t);
            R.dLambda_rad = dLambda;
            R.dLambda_deg = rad2deg(dLambda);
        end


        function dV = dv_combined(v1, v2, dI)
            % Vector addition formula from assignment
            dV = sqrt( v1^2 + v2^2 - 2*v1*v2*cos(dI) );
        end


        function lambda_s = interceptor_arrival_longitude(tw, eta, P)
            % Model used:
            % 1) coast in initial circular LEO for tw -> argument of latitude u increases
            % 2) Burn 1 occurs at that inertial position; inclination reduced by eta*i0
            % 3) Transfer ellipse: perigee at that position; propagate to apogee (f=pi)
            % 4) Compute inertial arrival position -> longitude atan2(y,x)

            % LEO argument of latitude after waiting (two-body, circular)
            u1 = P.u0 + P.nLEO*tw;

            % inclination after Burn 1 (since eta fraction of total i0 is removed)
            i1 = (1 - eta) * P.i0;

            % Set argument of perigee so that perigee direction aligns with current u1
            % For f=0 at perigee: u = omega + f = omega -> omega = u1
            omega = u1;

            % In perifocal frame, at apogee (f=pi): r = [-rA; 0; 0]
            r_pf = [-P.rGEO; 0; 0];

            % Transform perifocal -> inertial: R3(Omega)*R1(i1)*R3(omega)
            r_eci = Functioncollector.R3(P.Om) * Functioncollector.R1(i1) * Functioncollector.R3(omega) * r_pf;

            lambda_s = atan2(r_eci(2), r_eci(1));   % [-pi,pi]
        end


        function ang = wrapToPi(ang)
            % Wrap angle to (-pi, pi]
            ang = mod(ang + pi, 2*pi) - pi;
        end


        function R = R1(a)
            ca = cos(a); sa = sin(a);
            R = [ 1  0   0;
                  0  ca -sa;
                  0  sa  ca];
        end


        function R = R3(a)
            ca = cos(a); sa = sin(a);
            R = [ ca -sa  0;
                  sa  ca  0;
                  0   0   1];
        end

    end
end
