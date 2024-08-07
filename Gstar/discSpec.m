%
% Function: discSpec(par)
%
% Uses the continuous relaxation spectrum extracted using contSpec()
% to determine an approximate discrete approximation.
%
% Input: Communicated by the datastructure "par"
%
%        GstFile  = name of file that contains G*(w) in 3 columns [w Gp Gpp]
%                   default: 'Gst.dat' is assumed.
%        verbose  = 1, then prints onscreen messages, and prints datafiles
%        plotting = 1, then plots to stdio.
%
%
% Output: 
%         [g tau] = spectrum
%        
%         dmodes.dat : Prints the [g tau] for the particular Nopt
%         Nopt.dat   : Error and AIC for number of modes explored
%         Gfitd.dat  : The discrete G* for Nopt [w Gp Gpp] 
%

function [g tau] =  discSpec(par)

	%
    % Add appropriate subdirectories to search path
    %
	addpath('./disc:./common')

	if nargin == 0
		par = SetParameters();  % Load in global settings
	end

	if(par.verbose)
		fprintf('\n(*) Start\n(*) Loading Data Files: ...');
	end

	[w, Gexp, s, H, Nv, Gc, Cerror] = InitializeDiscSpec(par);

	n    = length(w);
	ns   = length(s);
	npts = length(Nv);

    % range of wtBaseDist scanned
    wtBase = par.deltaBaseWeightDist * [1: floor(1./par.deltaBaseWeightDist)];
    AICbst = zeros(length(wtBase), 1);
    Nbst   = zeros(length(wtBase), 1);  % nominal number of modes
    nzNbst = zeros(length(wtBase), 1);  % actual number of nonzero modes

	% main loop over wtBaseDist
	for ib = 1:length(wtBase)

		wb = wtBase(ib);
		
		% Find the distribution of nodes you need
		wt = GetWeights(H, w, s, wb);
		
		% Scan the range of number of Maxwell modes N = (Nmin, Nmax) 
		ev   = zeros(1, npts);
		nzNv = zeros(1, npts);  % number of nonzero modes 
		
		for i = 1:length(Nv)
			N = Nv(i);
			[z, hz] = GridDensity(log(s), wt, N);         % select "tau" Points
			[g, tau, ev(i), ~] = MaxwellModes(z, w, Gexp, par.plateau);    % get g*i
			nzNv(i) = length(g);
		end

		% store the best solution for this particular wb
		AIC           = 2 * Nv + 2 * Cerror * ev;
		AICbst(ib)    = min(AIC);
		[~, minIndex] = min(AIC);
		Nbst(ib)      = Nv(minIndex);
		nzNbst(ib)    = nzNv(minIndex);
	end

	% global best settings of wb and Nopt; note this is nominal Nopt (!= length(g) due to NNLS)
	[~, minIndex] = min(AICbst);
	Nopt          = round(Nbst(minIndex));
	wbopt         = wtBase(minIndex);

    %
    % Recompute the best data-set stats
    %
	wt = GetWeights(H, w, s, wbopt);
	[z, hz] = GridDensity(log(s), wt, Nopt); 							% Select "tau" Points
	[g, tau, err, cKp] = MaxwellModes(z, w, Gexp, par.plateau); % Get g_i, tau_i

	[succ, gf, tauf]   = FineTuneSolution(tau, w, Gexp, par.plateau);
	if succ
		g = gf;
		tau = tauf;
	end

	% Check if modes are close enough to merge
	[~, indx] = sort(tau);
	tau = tau(indx);
	tauSpacing = tau(2:end) ./ tau(1:end-1);
	itry = 0;

	if par.plateau
		g(1:end-1) = g(indx);
	else
		g = g(indx);
	end

	while min(tauSpacing) < par.minTauSpacing && itry < 3
		fprintf('\tTau Spacing < minTauSpacing\n');
		
		[~, imode] = min(tauSpacing); % merge modes imode and imode + 1
		[g, tau] = mergeModes_magic(g, tau, imode);
		[succ, g, tau] = FineTuneSolution(tau, w, Gexp, par.plateau);
		if succ
			g = gf;
			tau = tauf;
		end
		
		tauSpacing = tau(2:end) ./ tau(1:end-1);
		itry = itry + 1;
	end

	if par.plateau
		G0 = g(end);
		g  = g(1:end-1);
	end

	if par.verbose
		fprintf('(*) Number of optimum nodes = %d\n', length(g));
	end

	%
	% Some Plotting
	%

	if(par.plotting)

		subplot(2,1,1)

		loglog(tau, g,'o-')
		xlabel('tau')
		ylabel('g')
		title('discrete spectrum')

		subplot(2,1,2)
		if par.plateau
			PlotMaxwellModes(g, tau, w, Gexp(1:n), Gexp(n+1:end), G0)
		else
			PlotMaxwellModes(g, tau, w, Gexp(1:n), Gexp(n+1:end))
		end

		xlabel('w')
		ylabel('G*')
		hold off;

	end

	%
	% Some Printing
	%
	if(par.verbose)

		if par.plateau
			fprintf('(*) G0 = %e\n', G0);
		end

		f1 = fopen('output/dmodes.dat','w');
		fprintf('(*) Condition number of matrix equation: %e\n',cKp);

		fprintf('\n\t\tModes\n\t\t-----\n\n');
		fprintf('i \t    g(i) \t    tau(i)\n');
		fprintf('---------------------------------------\n');
		for i = 1:length(g)
			fprintf('%d \t %9.5e \t %9.5e\n',i, g(i), tau(i));
			fprintf(f1,'%d \t %9.5e \t %9.5e\n',i, g(i), tau(i));
		end
		fprintf('\n');
		fclose(f1);

		f2 = fopen('output/Nopt.dat','w');
		for i = 1:npts
			fprintf(f2,'%3d\t%e\t%e\n',Nv(i), ev(i), AIC(i));
		end
		fclose(f2);

		f3 = fopen('output/Gfitd.dat','w');
		K   = kernel(H, w, s);
		
		for i = 1:n
			fprintf(f3,'%e\t%e\t%e\n',w(i),K(i), K(n+i));
		end
		fclose(f3);

	end

end

% 
% Read Data Files and Initialize
%
function   [w, Gexp, s, H, Nv, Gc, Cerror] = InitializeDiscSpec(par)

	% load Exp Data
	[w, Gexp] = GetExpData(par.GstFile);
	n        = numel(w);
 
	% load CRS
    if par.plateau
    	data = dlmread('output/H.dat', '\t', 1, 0);
    else
    	data = load('output/H.dat');
    end
	s       = data(:,1);
	H       = data(:,2);
	ns      = numel(s);

	% range of N scanned
	Nmax = min(floor(3.0 * log10(max(w)/min(w))), n/4); % maximum N
	if par.MaxNumModes > 0
		Nmax = min(Nmax, par.MaxNumModes);
	end
	Nmin = max(floor(0.5 * log10(max(w)/min(w))), 3);   % minimum N
	Nv = Nmin:Nmax;

	% Estimate Error Weight from Continuous Curve Fit
	kernMat = getKernMat(s, w);


	if par.plateau
		try

			fid = fopen('output/H.dat', 'r');
			firstLine = textscan(fid, '%s', 'Delimiter', '\n');
			fclose(fid);

			firstLineTokens = strsplit(firstLine{1}{1});
			if length(firstLineTokens) > 2
				G0 = str2double(firstLineTokens{end});
			end
		catch
			error('Problem reading G0 from H.dat; Plateau = True');
		end
		Gc = kernel_prestore(H, kernMat, G0);
	else
		Gc = kernel_prestore(H, kernMat);
	end

	Cerror = 1 ./ (std((Gc ./ Gexp - 1)));  % Cerror = 1.?

end

