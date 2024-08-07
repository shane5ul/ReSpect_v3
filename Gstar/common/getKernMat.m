function kernMat = getKernMat(s, w)
    % Furnish kernMat which helps faster kernel evaluation, given s, w
    % Generates a 2n*ns matrix [(ws^2/1+ws^2) | (ws/1+ws)]' * hs, which can be 
    % multiplied with exp(H) to get predicted G*

    ns = length(s);
    hsv = zeros(ns, 1);

    hsv(1) = 0.5 * log(s(2) / s(1));
    hsv(ns) = 0.5 * log(s(ns) / s(ns-1));
    hsv(2:ns-1) = 0.5 * (log(s(3:ns)) - log(s(1:ns-2)));

    [S, W] = meshgrid(s, w);
    ws = S .* W;
    ws2 = ws .^ 2;

    kernMat = [ws2 ./ (1 + ws2); ws ./ (1 + ws2)] .* hsv';

end
