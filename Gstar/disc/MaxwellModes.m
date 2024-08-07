%
% Function: MaxwellModes(input)
%
% Solves the linear least squares problem to obtain the DRS
%
% Input: z = points distributed according to the density,
%        w = n*1 vector contains frequencies,
%        Gp, Gpp
%        Prune = Avoid modes with -ve weights (=1), or don't care (0) 
%
% Output: g, tau = spectrum 
%         error = relative error between the input data and the G*(w) inferred from the DRS
%         condKp = condition number
%
function [g, tau, err, condKp] = MaxwellModes(z, w, Gexp, isPlateau)

    N   = length(z);
    tau = exp(z);
    n   = length(w);
    
    % Prune small -ve weights g(i)
    [g, err, condKp] = nnLLS(w, tau, Gexp, isPlateau);
    
    % First remove runaway modes outside window with potentially large weight
    izero      = find(max(w)*min(tau) < 0.02 | min(w)*max(tau) > 50);
    tau(izero) = [];
    g(izero)   = [];
    
    % Search for small weights (gi) 
    if isPlateau
        izero = find(g(1:end-1)/max(g(1:end-1)) < 1e-8);
    else
        izero = find(g/max(g) < 1e-8);
    end
    tau(izero) = [];
    g(izero)   = [];
end

