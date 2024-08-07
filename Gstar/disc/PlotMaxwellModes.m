%
% Function PlotMaxwellModes(input)
%
% Plots experminetal dynamic moduli and dynamic moduli obtained from the
% corresponding CRS
%
% Input: g, t = spectrum
%        w = n*1 vector contains frequencies,
%        Gp, Gpp
%

function PlotMaxwellModes(g, t, w, Gp, Gpp, varargin);

	if nargin > 5
		G0 = varargin{1};
	else
		G0 = 0;
	end

	N  = length(g);

	[X, Y]  = meshgrid(t,w);
	ws      = X.*Y;
	ws2     = ws.^2;
	clear X Y;

	GpM  = (ws2./(1+ws2)) * g + G0;
	GppM = (ws./(1+ws2))  * g;

	loglog(w,Gp, 'bo', w, GpM, 'k-','LineWidth',2); hold on;
	loglog(w,Gpp,'bo', w, GppM,'k-','LineWidth',2);

end
