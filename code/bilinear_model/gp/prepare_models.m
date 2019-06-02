% Prepare the GP model used for control
disp('Loading GP model ...');

GPDATA = load(['results/' selected_gpmodel]);

GPMPCDATA = struct;
GPMPCDATA.X = GPDATA.X;
GPMPCDATA.Y = GPDATA.Y;
GPMPCDATA.hyp = GPDATA.hyp;
GPMPCDATA.cov = GPDATA.covfunc;
GPMPCDATA.lik = GPDATA.likfunc;
GPMPCDATA.input_model = GPDATA.input_model;
GPMPCDATA.normalization = GPDATA.normalization;

clear GPDATA;