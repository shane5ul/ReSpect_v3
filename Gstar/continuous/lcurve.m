function [lamM, lam, rho, eta, logP, Hlambda] = lcurve(Gexp, Hgs, kernMat, par, varargin)
    % Function: lcurve(input)
    %
    % Input: Gexp    = 2n*1 vector [G';G"],
    %        Hgs     = guessed H,
    %        kernMat = matrix for faster kernel evaluation
    %        par     = parameter dictionary
    %        G0      = optionally
    %
    % Output: lamC and 3 vectors of size npoints*1 contains a range of lambda, rho
    %         and eta. "Elbow"  = lamC is estimated using a *NEW* heuristic AND by Hansen method
    %
    % March 2019: starting from large lambda to small cuts calculation time by a lot
    %             also gives an error estimate 

    if(par.plateau)
        G0 = varargin{1};
    end

    npoints = round(par.lamDensity * (log10(par.lam_max) - log10(par.lam_min)));
    hlam    = (par.lam_max / par.lam_min)^(1 / (npoints - 1));
    lam     = par.lam_min * hlam .^ (0:(npoints-1));

    eta     = zeros(npoints, 1);
    rho     = zeros(npoints, 1);
    logP    = zeros(npoints, 1);

    H  = Hgs;
    n  = length(Gexp);
    ns = length(H);
    nl = ns - 2;

    logPmax = -inf; % so nothing surprises me!
    Hlambda = zeros(ns, npoints);

    % Error Analysis: Furnish A_matrix
    Amat    = getAmatrix(length(H));
    LogDetN = abs(logdet(Amat));

    % This is the costliest step
    for i = npoints:-1:1

        lamb = lam(i);

        if par.plateau
            [H, G0] = getH(lamb, Gexp, H, kernMat, G0);
            rho(i) = norm((1 - kernel_prestore(H, kernMat, G0) ./ Gexp));
            Bmat = getBmatrix(H, kernMat, Gexp, G0);
        else
            H = getH(lamb, Gexp, H, kernMat);
            rho(i) = norm((1 - kernel_prestore(H, kernMat) ./ Gexp));
            Bmat = getBmatrix(H, kernMat, Gexp);
        end


        eta(i) = norm(diff(H, 2));
        Hlambda(:, i) = H;

        LogDetC = abs(logdet(lamb * Amat + Bmat));
        V = rho(i)^2 + lamb * eta(i)^2;

        % This assumes a prior exp(-lam)
        logP(i) = -V + 0.5 * (LogDetN + ns * log(lamb) - LogDetC) - lamb;

        if logP(i) > logPmax
            logPmax = logP(i);
        elseif logP(i) < logPmax - 18
            break;
        end
    end

    % Truncate all to significant lambda
    lam     = lam(i:end);
    logP    = logP(i:end);
    eta     = eta(i:end);
    rho     = rho(i:end);
    logP    = logP - max(logP);
    Hlambda = Hlambda(:, i:end);

    % Currently using both schemes to get optimal lamC
    % New lamM works better with actual experimental data  
    plam = exp(logP);
    plam = plam / sum(plam);
    lamM = exp(log(lam)*plam);

    % Dialling in the Smoothness Factor
    if par.SmFacLam > 0
        lamM = exp(log(lamM) + par.SmFacLam * (max(log(lam)) - log(lamM)));
    elseif par.SmFacLam < 0
        lamM = exp(log(lamM) + par.SmFacLam * (log(lamM) - min(log(lam))));
    end
    
end

function A = getAmatrix(ns)
    % Generate symmetric matrix A = L' * L required for error analysis:
    % helper function for lcurve in error determination
    % L is a ns*ns tridiagonal matrix with 1, -2, and 1 on its diagonal

    nl = ns - 2;
    L = diag(ones(ns-1, 1), 1) + diag(ones(ns-1, 1), -1) + diag(-2 * ones(ns, 1));
    L = L(2:nl+1, :);
    
    A = L' * L;
end

function B = getBmatrix(H, kernMat, Gexp, varargin)
    % get the Bmatrix required for error analysis; helper for lcurve()
    %     not explicitly accounting for G0 in Jr because otherwise I get underflow problems

    n = floor(length(Gexp) / 2);  % Number of elements in Gexp (use floor for integer division)
    ns = size(H, 1);                % Number of rows in H
    nl = ns - 2;                     % Number of linear elements in H
    r = zeros(n, 1);                 % Initialize vector r

    % Relevant portion of Jacobian and residual
    Kmatrix = (1 ./ Gexp) .* ones(1, ns); 
    Jr = -kernelD(H, kernMat) .* Kmatrix;      % Jacobian Jr

    % Check for additional argument (G0)
    if nargin > 3
        G0 = varargin{1};
        r = (1 - kernel_prestore(H, kernMat, G0) ./ Gexp);
    else
        r = (1 - kernel_prestore(H, kernMat) ./ Gexp);
    end

    % B matrix calculation
    B = Jr' * Jr + diag(r' * Jr);  % Use diag for diagonal matrix

    return;
end