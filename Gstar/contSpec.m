% 
% Octave Note: Install and load 'optim' package:
% pkg install optim-1.6.2.tar.gz   % to install
% pkg load optim                   % to load in current session
%
% see https://github.com/gnu-octave/pkg?tab=readme-ov-file for more details
% on package management
%
%

% Function: contSpec(par)
%
% Using a simplified L-curve method to compute the continuous relaxation 
% spectra H(s) given G*(w) from an input file
%
% Uses parameters from SetParameters.m file:
% These parameters are prefixed by par.XXXXX, where XXXXX is the parameter
% Definitions of these parameters are provided in the SetParameters.m file
%
% Output:   Yields the fitted H(s), and optionally critical lambda "lamC"
%           If verbose is on: rho-eta.dat (used to determine lamC)
%                             H.dat       (the spectrum)
%                             Gfit.dat    (G* implied by the extracted H(s))
%
function  [H, lamC] = contSpec(par)

    %
    % Add appropriate subdirectories to search path
    %
	addpath('./continuous:./common');
	
	if nargin == 0
		par = SetParameters();  % Load in global settings
	end

	%
	% Load the data in
	% 
	if(par.verbose)
		fprintf('\n(*) Start\n(*) Loading Data File: %s...',par.GstFile);
	end

	[w, Gexp] = GetExpData(par.GstFile);

	if(par.verbose)
		fprintf('done\n(*) Initial Set up...');
	end

	tic

	%
	% Set up some internal variables
	%
	n    = length(w);
	ns   = par.ns;    % discretization of 'tau'

	wmin = w(1);
	wmax = w(n);

	switch par.FreqEnd
	case 1
		smin = exp(-pi/2)/wmax; smax = exp(pi/2)/wmin;
	case 2
		smin = 1/wmax; smax = 1/wmin;
	case 3
		smin = exp(+pi/2)/wmax; smax = exp(-pi/2)/wmin;
	end 

	hs   = (smax/smin)^(1/(ns-1));
	s    = smin * hs.^[0:ns-1]';

	kernMat = getKernMat(s, w);

    % get an initial guess for Hgs, G0
    if par.plateau
        [Hgs, G0]  = InitializeH(Gexp, w, s, kernMat, min(Gexp));
    else
        Hgs        = InitializeH(Gexp, w, s, kernMat);
	end

	t = toc;
	tic

	%
	% Find Optimum Lambda with 'lcurve'
	%
	if(par.lamC == 0)

	    if par.plateau
            [lamC, lam, rho, eta, logP, Hlam] = lcurve(Gexp, Hgs, kernMat, par, G0);
        else
            [lamC, lam, rho, eta, logP, Hlam] = lcurve(Gexp, Hgs, kernMat, par);
		end

	else
		lamC = par.lamC;
	end

	t = toc;

	if(par.verbose)
		fprintf('lamC = %12.4e (%5.1f seconds)\n(*) Extracting the continuous spectrum, ...',lamC, t);
	end

	tic;
	
	% Get the spectrum
    if par.plateau
        [H, G0] = getH(lamC, Gexp, Hgs, kernMat, G0);
        fprintf('G0 = %8.4e ...',G0)
    else
        H = getH(lamC, Gexp, Hgs, kernMat);
    end
	
	
	t = toc;

	%
	% Print some datafiles
	%
	if(par.verbose)
		fprintf('done (%5.1f seconds)\n(*) Writing and Printing, ...',t);
	end

	if(par.verbose)

		if(par.lamC == 0)
			f1  = fopen('output/rho-eta.dat','w');
			for i = 1:length(lam)
				fprintf(f1,'%e\t%e\t%e\n',lam(i),rho(i),eta(i));
			end
			fclose(f1);

			f1  = fopen('output/logP.dat','w');
			for i = 1:length(lam)
				fprintf(f1,'%e\t%e\n',lam(i),logP(i));
			end
			fclose(f1);

		end

		f2  = fopen('output/H.dat','w');

		if par.plateau
			fprintf(f2, '# G0 = %e\n', G0);
			K   = kernel_prestore(H, kernMat, G0);
		else
			K   = kernel_prestore(H, kernMat);
        end

		for i = 1:ns
			fprintf(f2,'%e\t%e\n', s(i), H(i));
		end
		fclose(f2);

		f3  = fopen('output/Gfit.dat','w');
		for i = 1:n
			fprintf(f3,'%e\t%e\t%e\n',w(i),K(i),K(n+i));
		end
		fclose(f3);


	end

	%
	% Graphing
	%    

	if(par.plotting)

		subplot(2,1,1)
		semilogx(s,H,'o-')
		xlabel('s')
		ylabel('H(s)')

		subplot(2,1,2);

		if par.verbose == 0
			if par.plateau
				K   = kernel_prestore(H, kernMat, G0);
			else
				K   = kernel_prestore(H, kernMat);
			end
		end

		loglog(w,Gexp(1:n),'o',w,K(1:n),'k-',w,Gexp(n+1:2*n),'s',w,K(n+1:2*n),'k-');
		xlabel('w')
		ylabel('G*(exp), G*(fit)')

	end

	if(par.verbose)
		fprintf('done\n(*) End\n');
	end

end
