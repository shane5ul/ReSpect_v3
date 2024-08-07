function [success, g, tau] = FineTuneSolution(tau, w, Gexp, isPlateau)
    % Given a spacing of modes tau, tries to do NLLS to fine tune it further
    % If it fails, then it returns the old tau back	   
    % Uses helper function: res_wG which computes residuals

    success = false;
    initError = norm(res_wG(tau, w, Gexp, isPlateau));
    

    try
        options = optimset('MaxIter', 5000, 'Display', 'off');
        lb = 0.02 / max(w) * ones(size(tau));
        ub = 50 / min(w) * ones(size(tau));
        [tau, ~, ~, exitflag] = lsqnonlin(@(x) res_wG(x, w, Gexp, isPlateau), tau, lb, ub, options);
        
        if exitflag > 0
            success = true;
        end
    catch
        % If optimization fails, keep the original tau
    end
    
    [g, tau, ~, ~] = MaxwellModes(log(tau), w, Gexp, isPlateau);   % Get g_i, taui
    finalError = norm(res_wG(tau, w, Gexp, isPlateau));
    
    % keep fine tuned solution, only if it improves things
    if finalError > initError
        success = false;
    end
end

% Helper function for final optimization problem
function residual = res_wG(tau, wexp, Gexp, isPlateau)
    [g, ~, ~]     = nnLLS(wexp, tau, Gexp, isPlateau);
    Gmodel        = zeros(size(Gexp));
    
    [S, W] = meshgrid(tau, wexp);
    ws = S .* W;
    ws2 = ws.^2;
    K = [ws2 ./ (1 + ws2); ws ./ (1 + ws2)];   % 2n * nmodes
    
    % add G0
    if isPlateau
        Gmodel = K * g(1:end-1);
        n = length(Gexp) / 2;
        Gmodel(1:n) = Gmodel(1:n) + g(end);	
    else
        Gmodel = K * g;
    end
    
    residual = (Gmodel ./ Gexp - 1);
end
