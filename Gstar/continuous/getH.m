function [H, G0] = getH(lam, Gexp, H, kernMat, varargin)

    % minimize_H  V(lambda) := ||(Gexp - kernel(H)) * (wexp/Gexp)||^2 +  lambda * ||L H||^2

    % Input : lambda  = regularization parameter ,
    %         Gexp    = experimental data,
    %         H       = guessed H,
    %         kernMat = matrix for faster kernel evaluation
    %         G0      = optional

    % Output : H_lam, [G0]
    %         Default uses Trust-Region Method with Jacobian supplied by jacobianLM
    %

    % Check if there is an additional argument (G0)
    options = optimset('Jacobian', 'on', 'Display', 'off');
    if (nargout == 2)
        Hplus   = [H; varargin{1}];
        res_lsq = lsqnonlin(@(x) getResidualJacobian(x, lam, Gexp, kernMat), Hplus, [], [], options);
        H       = res_lsq(1:end-1);
        G0      = res_lsq(end);
    else
        res_lsq = lsqnonlin(@(x) getResidualJacobian(x, lam, Gexp, kernMat), H, [], [], options);
        H       = res_lsq;
    end

end

%
% HELPER FUNCTION: Gets Residuals r and Jacobian J
% If called with only one output variable, only returns residual
%
function [r Jr] = getResidualJacobian(H, lam, Gexp, kernMat)

    n  = floor(size(kernMat, 1) / 2);
    ns = size(kernMat, 2);
	nl  = ns - 2;

	% Get the residual vector first
	% r = vector of size (2n+nl,1) regardless
	r   = zeros(2*n + nl,1);

	% Check for G0
    if length(H) > ns
        plateau  = 1;
        G0       = H(end);
        H        = H(1:end-1);
        r(1:2*n) = (1 - kernel_prestore(H, kernMat, G0) ./ Gexp);
    else
        plateau  = 0;
        r(1:2*n) = (1 - kernel_prestore(H, kernMat) ./ Gexp);
    end

    % Curvature constraint
    r(2*n+1:end) = sqrt(lam) * diff(H, 2);  % Second derivative

	% Furnish the Jacobian Jr
	% size depends on plateau on or off
	if(nargout == 2)

        % L is a nl*ns tridiagonal matrix with 1 -2 and 1 on its diagonal.
        L = diag(ones(1, ns-1), 1) + diag(ones(1, ns-1), -1) + diag(-2 * ones(1, ns));
        L = L(2:nl+1, :);
        
        % Furnish the Jacobian Jr - (2n+nl)*ns matrix
        % Kmatrix is 2*n * ns matrix
        Kmatrix = (1 ./ Gexp) * ones(1, ns);

        if plateau
            Jr = zeros(2*n + nl, ns+1);
            
            Jr(1:2*n, 1:ns) = -kernelD(H, kernMat) .* Kmatrix;
            Jr(1:n, ns+1) = -1 ./ Gexp(1:n);  % nonzero dr_i/dG0 only for G'
            Jr(2*n+1:2*n+nl, 1:ns) = sqrt(lam) * L;
            Jr(2*n+1:2*n+nl, ns+1) = zeros(nl, 1);  % column for dr_i/dG0 = 0
        else
            Jr = zeros(2*n + nl, ns);
            Jr(1:2*n, 1:ns) = -kernelD(H, kernMat) .* Kmatrix;
            Jr(2*n+1:2*n+nl, 1:ns) = sqrt(lam) * L;
        end

	end


end

