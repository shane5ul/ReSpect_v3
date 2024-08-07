%  Function: kernelD(input)
%
%  outputs the 2n*ns dimensional vector DK(H)(w)
%  It approximates dK_i/dH_j = K * e(H_j):
%
%  Input: H       = substituted CRS,
%         kernMat = matrix for faster kernel evaluation
%
%  Output: DK     = Jacobian of H

function DK = kernelD(H, kernMat)
    [n, ns] = size(kernMat);
    n = n / 2;
    
    Hsuper = ones(2*n, 1) * exp(H(:)');
    DK = kernMat .* Hsuper;
end
