
%
% Set parameters for calculation
%
% August 5, 2024: Using LLMs to port pyReSpect-freq algorithm
%                 Simplifying structure greatly
%


function par = SetParameters()

				%=====================================%
				% O V E R A L L   P A R A M E T E R S %
				%=====================================%

	%
	% Default filename for Gp and Gpp data
	% It should contain G*(w) in 3 columns [w Gp Gpp]
	%
	par.GstFile  = 'Gst.dat';

	%
	% Plateau Modulus: is there a residual plateau in the data? ON (=1) or OFF (0)
	%
	par.plateau  = 0;

	%
	% Printing to screen and files is ON (=1) or OFF (0)
	%
	par.verbose  = 1;

	%
	% Plotting functions are ON (=1) or OFF (0)
	%
	par.plotting = 1;

				%=======================================%
				% C O N T I N U O U S   S P E C T R U M %
				%=======================================%
	%
	% Number of grid points to represent the continuous spectrum
	%
	par.ns = 100;

	%
	% Treatment of frequency window ends:
	%
	%  = 3 : t = 1/w - strict condition
	%  = 2 : t = 1/w
	%  = 1 : t = 1/w + lenient condition
	%
	par.FreqEnd = 1;


	%
	% Specify lambda_C instead of using the one inferred from the L-curve
	% If 0, then use Bayesian inference to determine lambda.
	%
	par.lamC = 0;

	%
	% Smoothing Factor: Indirect way of controlling lambda_C, relative to
	% the one inferred from L-curve.
	%
	% Set between -1 (lowest lambda explored) and 1 (highest lambda explored);
	% When set to 0, using lambda_C determined from the L-curve
	%
	par.SmFacLam = 0;

	%
	% Max Number of Modes: Set = 0, if you want it to automatically determine
	% the optimal number of modes; otherwise this parameter sets the maximum possible.
	par.MaxNumModes  = 0; 

	%
	% contSpec: limits of lambda explored
	%           
	par.lam_min = 1e-10;
	par.lam_max = 1e3;

	%
	% Lambda Density per Decade for Lcurve: set = 2 or more
	%
	par.lamDensity  = 3;

	%
	% discSpec: how finely to sample BaseWeightDist
	%           smaller value increases cost; typical range 0.10 - 0.25
	par.deltaBaseWeightDist = 0.1;

	% discSpec: how close do successive modes (tau2/tau1) have to be before we 
	%           try to mege them
	par.minTauSpacing       = 1.25;

end

