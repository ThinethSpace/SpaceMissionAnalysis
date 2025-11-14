classdef Orbit
    % Lightweight orbital utilities
    methods (Static)
        function [t_vec, R_eci, V_eci] = propagateKeplerNewton(a,e,i,Omega,omega,nu0,mu,t0,tf,dt)
            % Unperturbed Keplerian propagation using Newton's method.
            % Angles in radians. Returns ECI position/velocity at each time step.

            % nu0 -> E0 -> M0
            E0 = 2*atan2( sqrt(1-e)*sin(nu0/2), sqrt(1+e)*cos(nu0/2) );
            M0 = E0 - e*sin(E0);

            n  = sqrt(mu/a^3);                     % mean motion
            Nt = floor((tf - t0)/dt) + 1;
            t_vec = linspace(t0, tf, Nt).';
            R_eci = zeros(Nt,3);
            V_eci = zeros(Nt,3);

            % Fixed rotation PQW->ECI (Keplerian = elements constant)
            cO = cos(Omega); sO = sin(Omega);
            ci = cos(i);     si = sin(i);
            co = cos(omega); so = sin(omega);
            Q = [ cO*co - sO*so*ci,   -cO*so - sO*co*ci,   sO*si;
                  sO*co + cO*so*ci,   -sO*so + cO*co*ci,  -cO*si;
                  so*si,               co*si,             ci     ];

            p = a*(1 - e^2);

            for k = 1:Nt
                t = t_vec(k);
                M = M0 + n*(t - t0);               % mean anomaly

                % Newton for Kepler's equation: E - e sinE = M
                E = Orbit.wrapToPiLocal(M);        % good initial guess
                for it = 1:10
                    f  = E - e*sin(E) - M;
                    fp = 1 - e*cos(E);
                    dE = -f/fp;
                    E  = E + dE;
                    if abs(dE) < 1e-12, break; end
                end

                % True anomaly, radius
                nu = 2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );
                r  = a*(1 - e*cos(E));

                % PQW state
                r_pqw = [ r*cos(nu); r*sin(nu); 0 ];
                v_pqw = sqrt(mu/p) * [ -sin(nu); e + cos(nu); 0 ];

                % To ECI
                R_eci(k,:) = (Q * r_pqw).';
                V_eci(k,:) = (Q * v_pqw).';
            end
        end
    end

    methods (Static)
        function [t_vec, R_eci, V_eci, Omegadot] = propagateKeplerJ2Node( ...
                a,e,i,Omega0,omega,nu0,mu,J2,Re,t0,tf,dt)
            n  = sqrt(mu/a^3);
            p  = a*(1 - e^2);
            Omegadot = -1.5 * J2 * (Re^2 / p^2) * n * cos(i);  % J2 RAAN drift
    
            E0 = 2*atan2( sqrt(1-e)*sin(nu0/2), sqrt(1+e)*cos(nu0/2) );
            M0 = E0 - e*sin(E0);
    
            Nt = floor((tf - t0)/dt) + 1;
            t_vec = linspace(t0, tf, Nt).';
            R_eci = zeros(Nt,3); V_eci = zeros(Nt,3);
    
            for k = 1:Nt
                t = t_vec(k);
                Omega = Omega0 + Omegadot*(t - t0);   % Ω(t)
    
                M = M0 + n*(t - t0);
                E = mod(M + pi, 2*pi) - pi;           % initial guess
                for it = 1:10
                    dE = -(E - e*sin(E) - M) / (1 - e*cos(E));
                    E = E + dE;
                    if abs(dE) < 1e-12, break; end
                end
    
                nu = 2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );
                r  = a*(1 - e*cos(E));
                r_pqw = [ r*cos(nu); r*sin(nu); 0 ];
                v_pqw = sqrt(mu/p) * [ -sin(nu); e + cos(nu); 0 ];
    
                cO = cos(Omega); sO = sin(Omega);
                ci = cos(i);     si = sin(i);
                co = cos(omega); so = sin(omega);
                Q = [ cO*co - sO*so*ci,   -cO*so - sO*co*ci,   sO*si;
                      sO*co + cO*so*ci,   -sO*so + cO*co*ci,  -cO*si;
                      so*si,               co*si,             ci     ];
    
                R_eci(k,:) = (Q*r_pqw).';
                V_eci(k,:) = (Q*v_pqw).';
            end
        end
    
        function i = ssoInclination(a,e,OmegaDot_desired,mu,J2,Re)
            n  = sqrt(mu/a^3);  p = a*(1 - e^2);
            cosi = -OmegaDot_desired / (1.5 * J2 * (Re^2/p^2) * n);
            cosi = max(-1,min(1,cosi));
            i = acos(cosi);
        end
    end
    
        methods (Static, Access = private)
            function ang = wrapToPiLocal(theta)
                % Wrap angle(s) to (-pi, pi]
                ang = mod(theta + pi, 2*pi) - pi;
            end
        end
        methods (Static)
            function [t_vec, R_eci, V_eci] = propagateRK4(t0, tf, dt, r0, v0, accFcn)
                Nt    = floor((tf - t0)/dt) + 1;
                t_vec = linspace(t0, tf, Nt).';
                R_eci = zeros(Nt, 3);
                V_eci = zeros(Nt, 3);
    
                y = [r0(:); v0(:)];   % 6x1
                R_eci(1,:) = y(1:3).';
                V_eci(1,:) = y(4:6).';
    
                for k = 1:Nt-1
                    t = t_vec(k);
                    h = dt;
    
                    r = y(1:3); v = y(4:6);
                    a1 = accFcn(t, r, v);
                    k1 = [v; a1];
    
                    r = y(1:3) + 0.5*h*k1(1:3); v = y(4:6) + 0.5*h*k1(4:6);
                    a2 = accFcn(t + 0.5*h, r, v);
                    k2 = [v; a2];
    
                    r = y(1:3) + 0.5*h*k2(1:3); v = y(4:6) + 0.5*h*k2(4:6);
                    a3 = accFcn(t + 0.5*h, r, v);
                    k3 = [v; a3];
    
                    r = y(1:3) + h*k3(1:3); v = y(4:6) + h*k3(4:6);
                    a4 = accFcn(t + h, r, v);
                    k4 = [v; a4];
    
                    y = y + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    
                    R_eci(k+1,:) = y(1:3).';
                    V_eci(k+1,:) = y(4:6).';
                end
            end
        end
end
