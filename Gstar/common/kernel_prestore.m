function result = kernel_prestore(H, kernMat, varargin)

    % 08/05/2024: Using chatGPT to convert!

    %  turbocharging kernel function evaluation by prestoring kernel matrix
    %  Date    : 8/17/2018
    %  Function: kernel_prestore(input) returns K*h, where h = exp(H)
        
    %  Same as kernel, except prestoring hs, S, and W to improve speed 3x.
    
    %  outputs the 2n*1 dimensional vector K(H)(w) which is comparable to G* = [G'|G"]'
    %  3/11/2019: returning Kh + G0
            
    %  Input: H = substituted CRS,
    %         kernMat = 2n*ns matrix [(ws^2/1+ws^2) | (ws/1+ws)]'*hs

    if nargin > 2
        n = floor(size(kernMat, 1) / 2);
        G0v = zeros(2 * n, 1);
        G0v(1:n) = varargin{1};
    else
        G0v = 0;
    end

    result = kernMat * exp(H) + G0v;
end
