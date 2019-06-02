clear all
close all

%% Load data

dist_type = 'est';

[xTrain, yTrain] = gp_load_data('train', dist_type);

use_normalization = true;
signals = gp_create_signals(xTrain, yTrain, use_normalization);


%% Prepare data for GP model
signals_model = SignalsModel(signals);

% Select features for training
% Output is 'temp', which can have auto-regression
%{
input_model = {... %     'tod', ... %     'dow', ... 
    {'temp', [2 1]}, ...
    {'d1', [2 1 0]}, ...
    {'d2', [2 1 0]}, ...
    {'d3', [2 1 0]}, ...  % {'d4', 4:-1:0}, ... 
    {'d5', [2 1 0]}, ...
    {'d6', [2 1 0]}, ...  
    {'d7', [2 1]}, ...
    {'blind', [2 1 0]}, ...
    {'light', [2 1 0]}, ... 
    {'cooling', [2 1 0]}, ... 
    {'heating', [2 1 0]}, ...
    };
%}
input_model = {... %     'tod', ... %     'dow', ... 
    {'temp', 4:-1:1}, ...
    {'d1', 4:-1:0}, ...
    {'d2', 4:-1:0}, ...
    {'d3', 4:-1:0}, ...  % {'d4', 4:-1:0}, ... 
    {'d5', 4:-1:0}, ...
    {'d6', 4:-1:0}, ...  
    {'d7', 4:-1:0}, ...
    {'blind', 4:-1:0}, ...
    {'light', 4:-1:0}, ... 
    {'cooling', 4:-1:0}, ... 
    {'heating', 4:-1:0}, ...
    };

signals_model.setInputs(input_model);
signals_model.setOutput('temp');

% Obtain all IO data
[X_all, Y_all] = signals_model.getIOVectors([], 'removeNaNs', true);

% Select indices for training
N = size(Y_all,1);
trainIdx = 1:5000;  % (N-3000):N;
Ntrain = length(trainIdx);
X_train = X_all(trainIdx,:);
Y_train = Y_all(trainIdx,:);

% Clean up memory a bit
clear xTrain yTrain


%% Train GP model

% Dimension of features
D = size(X_train, 2);

% Select the covariance function, mean function, etc.
hyp = struct;
covfunc = 'covSEard';  hyp.cov = zeros(D+1,1);  covfuncshort = 'se';
likfunc = 'likGauss';  hyp.lik = 0;
meanfunc = [];  hyp.mean = [];  meanfuncshort = 'meanzero';

% Use minFunc for faster optimization; if not available, use minimize()
% hyp = minimize_minfunc(hyp, @gp, -1000, @infExact, meanfunc, covfunc, likfunc, X_train, Y_train);

hyp = minimize(hyp, @gp, 2000, @infExact, meanfunc, covfunc, likfunc, X_train, Y_train);


%% Only for ARD kernel, print the names of features that have large
% lengthscales
% disp('Features with large lengthscales are:');
% maxhyp = 4.3;
% irrelevant_features = what_is_the_input(input_model, find(hyp.cov(1:end-1) > maxhyp))


%% Save model
RESULTS = struct;
RESULTS.covfunc = covfunc;
RESULTS.likfunc = likfunc;
RESULTS.meanfunc = meanfunc;
RESULTS.covfuncshort = covfuncshort;
RESULTS.meanfuncshort = meanfuncshort;
RESULTS.hyp = hyp;
RESULTS.Ntrain = Ntrain;
RESULTS.input_model = input_model;
RESULTS.normalization = signals.getSignalNormalization();
RESULTS.D = D;
RESULTS.X = X_train;
RESULTS.Y = Y_train;
RESULTS.dist_type = dist_type;

matname = sprintf('results/gp-%s-%d-%s-%s-%s', dist_type, Ntrain, covfuncshort, meanfuncshort, 'lag4'); %datestr(now(), 'yyyymmdd-HHMMSS'));
save(matname, '-struct', 'RESULTS');


%% Validation
performance = test_performance(matname, RESULTS.normalization, true);
