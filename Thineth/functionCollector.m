classdef functionCollector
    % All helper functions as static methods in a single file.

    methods(Static)

        %% ===== Horizons reader =====
        function S = readHorizonsCartesian(filepath)
            % Reads JPL Horizons cartesian states between $$SOE and $$EOE.
            % Output:
            %   S.jd   [Nx1] Julian date (days)
            %   S.date [Nx1] string calendar date
            %   S.r    [Nx3] AU
            %   S.v    [Nx3] AU/day

            txt = fileread(filepath);
            lines = splitlines(string(txt));

            iSOE = find(contains(lines,"$$SOE"),1,'first');
            iEOE = find(contains(lines,"$$EOE"),1,'first');
            if isempty(iSOE) || isempty(iEOE) || iEOE <= iSOE
                error("Could not find $$SOE/$$EOE in file: %s", filepath);
            end

            dataLines = lines(iSOE+1:iEOE-1);
            mask = contains(dataLines,",");
            dataLines = dataLines(mask);

            jd   = [];
            date = strings(0,1);
            r    = [];
            v    = [];

            for k = 1:numel(dataLines)
                L = strtrim(dataLines(k));
                if L == "" || startsWith(L,"***")
                    continue;
                end
                parts = split(L,",");

                if numel(parts) < 8
                    continue;
                end

                jd_k = str2double(strtrim(parts(1)));
                if isnan(jd_k)
                    continue;
                end
                date_k = strtrim(parts(2));

                x  = str2double(strtrim(parts(3)));
                y  = str2double(strtrim(parts(4)));
                z  = str2double(strtrim(parts(5)));
                vx = str2double(strtrim(parts(6)));
                vy = str2double(strtrim(parts(7)));
                vz = str2double(strtrim(parts(8)));

                if any(isnan([x y z vx vy vz]))
                    continue;
                end

                jd(end+1,1) = jd_k; %#ok<AGROW>
                date(end+1,1) = date_k; %#ok<AGROW>
                r(end+1,:) = [x y z]; %#ok<AGROW>
                v(end+1,:) = [vx vy vz]; %#ok<AGROW>
            end

            S = struct('jd',jd,'date',date,'r',r,'v',v);
        end

        %% ===== Lambert solver (secant in a) =====
        function [v1,v2,a,info] = lambertSecant_a(r1,r2,dt,mu,opts)
            arguments
                r1 (3,1) double
                r2 (3,1) double
                dt (1,1) double {mustBePositive}
                mu (1,1) double {mustBePositive}
                opts.longway (1,1) logical = false
                opts.maxIter (1,1) double {mustBeInteger,mustBePositive} = 60
                opts.tolRel (1,1) double {mustBePositive} = 1e-10
                opts.verbose (1,1) logical = false
            end

            r1m = norm(r1); r2m = norm(r2);
            cosd = dot(r1,r2)/(r1m*r2m);
            cosd = max(-1,min(1,cosd));
            dtheta = acos(cosd);
            if opts.longway
                dtheta = 2*pi - dtheta;
            end

            c = norm(r2 - r1);
            s = 0.5*(r1m + r2m + c);
            am = s/2;

            tof = @(a_) functionCollector.tofEllipse_fromA(a_,s,c,mu,dtheta);

            a1 = am * 1.0001;
            a2 = am * 1.2;

            f1 = tof(a1) - dt;
            f2 = tof(a2) - dt;

            expandCount = 0;
            while (~isfinite(f1) || ~isfinite(f2) || sign(f1)==sign(f2)) && expandCount < 30
                a2 = a2 * 1.4;
                f2 = tof(a2) - dt;
                expandCount = expandCount + 1;
            end

            converged = false;
            a = NaN;

            for k = 1:opts.maxIter
                denom = (f2 - f1);
                if abs(denom) < eps
                    break;
                end

                a3 = a2 - f2*(a2 - a1)/denom;
                if a3 <= am
                    a3 = 0.5*(a2 + am);
                end

                f3 = tof(a3) - dt;

                if opts.verbose
                    fprintf("iter %02d: a=%g, f=%g\n",k,a3,f3);
                end

                if abs(f3)/dt < opts.tolRel
                    a = a3;
                    converged = true;
                    break;
                end

                a1 = a2; f1 = f2;
                a2 = a3; f2 = f3;
            end

            if ~converged
                v1 = nan(3,1); v2 = nan(3,1);
                info = struct('converged',false,'dtheta',dtheta,'c',c,'s',s,'am',am);
                return;
            end

            [alpha,beta] = functionCollector.alphaBetaEllipse(a,s,c);

            p = (4*a*(s - r1m)*(s - r2m) / (c^2)) * (sin(0.5*(alpha + beta))^2);

            f    = 1 - (r2m/p)*(1 - cos(dtheta));
            g    = (r1m*r2m*sin(dtheta)) / sqrt(mu*p);
            gdot = 1 - (r1m/p)*(1 - cos(dtheta));

            if abs(g) < 1e-14
                error("Lambert: g is ~0 (degenerate geometry).");
            end

            v1 = (r2 - f*r1)/g;
            v2 = (gdot*r2 - r1)/g;

            info = struct('converged',true,'dtheta',dtheta,'a',a,'p',p,'iters',k);
        end

        function t = tofEllipse_fromA(a,s,c,mu,dtheta) %#ok<INUSD>
            am = s/2;
            if a <= am
                t = NaN; return;
            end
            [alpha,beta] = functionCollector.alphaBetaEllipse(a,s,c);
            t = sqrt(a^3/mu) * ( (alpha - sin(alpha)) - (beta - sin(beta)) );
        end

        function [alpha,beta] = alphaBetaEllipse(a,s,c)
            arg1 = s/(2*a);
            arg2 = (s - c)/(2*a);
            arg1 = max(0,min(1,arg1));
            arg2 = max(0,min(1,arg2));
            alpha = 2*asin(sqrt(arg1));
            beta  = 2*asin(sqrt(arg2));
        end

    end
end
