% 
% Helper subfunction which does the actual LLS problem
% helps MaxwellModes
%
function [g, err, condKp] = nnLLS(w, tau, Gexp, isPlateau)
    
    n      = length(Gexp) / 2;
    ntau   = length(tau);
    [S, W] = meshgrid(tau, w);

    ws  = S .* W;
    ws2 = ws .^ 2;
    K   = [ws2 ./ (1 + ws2); ws ./ (1 + ws2)]; % 2n * nmodes

    % K is n*ns [or ns+1]
    if isPlateau
        K = [K, ones(length(Gexp), 1)]; % G' needs some G0
        K(n+1:end, ntau+1) = 0;         % G" doesn't have G0 contrib
    end

    % gets (Gst/GstE - 1)^2, instead of  (Gst -  GstE)^2
    Kp     = diag(1 ./ Gexp) * K;
    condKp = cond(Kp);
    g      = lsqnonneg(Kp, ones(size(Gexp))); % MATLAB's lsqnonneg is similar to nnls in Python
    
    GstM  = K * g;
    err   = sum(((GstM ./ Gexp - 1)).^2);
end
