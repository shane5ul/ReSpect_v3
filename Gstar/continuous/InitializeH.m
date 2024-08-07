%
% Function: InitializeH(input)
%  Input:  Gexp    = 2n*1 vector [G';G"],
%          w       = experimental frequencies,
%          s       = relaxation modes,
%          kernMat = matrix for faster kernel evaluation
%          G0      = optional; if plateau is nonzero	

%  Output: H = guessed H
%           G0 = optional guess if *argv is nonempty	

function [Hlam, G0] = InitializeH(Gexp, w, s, kernMat, varargin)

    % To guess spectrum, pick a negative Hgs and a large value of lambda to get a
    % solution that is most determined by the regularization
    % March 2019; a single guess is good enough now, because going from large lambda to small
    %             lambda in lcurve.

	H       = -5.0 + sin(pi * s);
	lambda  = 1e0;

    if nargout > 1
        G0 = varargin{1};
        [Hlam, G0] = getH(lambda, Gexp, H, kernMat, G0);
    else
        Hlam = getH(lambda, Gexp, H, kernMat);
    end

end
