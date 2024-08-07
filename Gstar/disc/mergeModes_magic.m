function [newg, newtau] = mergeModes_magic(g, tau, imode)
    % Merge modes imode and imode+1 into a single mode
    % Return newg and newtau corresponding to this new mode

    iniGuess = [g(imode) + g(imode+1), 0.5 * (tau(imode) + tau(imode+1))];
    options = optimset('Display', 'off');
    res = fminsearch(@(par) costFcn_magic(par, g, tau, imode), iniGuess, options);

    newtau = tau;
    newtau(imode+1) = [];
    newtau(imode) = res(2);

    newg = g;
    newg(imode+1) = [];
    newg(imode) = res(1);
end

% Helper function for mergeModes; establishes cost function to minimize
function cost = costFcn_magic(par, g, tau, imode)
    gn = par(1);
    taun = par(2);

    g1 = g(imode);
    g2 = g(imode+1);
    tau1 = tau(imode);
    tau2 = tau(imode+1);

    wmin = min(1./tau1, 1./tau2) / 10;
    wmax = max(1./tau1, 1./tau2) * 10;

    % using clenshaw curtis for octave; may have to use something else in matlab
    cost = quadcc(@(w) normKern_magic(w, gn, taun, g1, tau1, g2, tau2), wmin, wmax);
end

% Helper function for costFcn and mergeModes
function normK = normKern_magic(w, gn, taun, g1, tau1, g2, tau2)
    wt = w * taun;
    Gnp = gn * (wt.^2) ./ (1 + wt.^2);
    Gnpp = gn * wt ./ (1 + wt.^2);

    wt = w * tau1;
    Gop = g1 * (wt.^2) ./ (1 + wt.^2);
    Gopp = g1 * wt ./ (1 + wt.^2);

    wt = w * tau2;
    Gop = Gop + g2 * (wt.^2) ./ (1 + wt.^2);
    Gopp = Gopp + g2 * wt ./ (1 + wt.^2);

    normK = (Gnp ./ Gop - 1).^2 + (Gnpp ./ Gopp - 1).^2;
end