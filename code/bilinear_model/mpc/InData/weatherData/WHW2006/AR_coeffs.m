% AR_coffs.m
% 
% stores AR coeffs

%% WHW 2006
%  AR of order 4 

% RGS raw data, predhor = 24
-0.045264917827551
   0.023166328048953
  -0.068075467934311
   0.685448697438124
  -4.700563969468201
   
sigma = 77.909579020994158

ACF_RGS(2) = -5.117865551199228e-04;
bounds = 0.004376923020075;

% TA raw data, predhor = 24
-0.006658251730190
  -0.062359836733178
  -0.088413137437204
   1.037095643071405
   0.168355788215680

sigma = 0.996752980816940
ACF_TA(2) = 2.340666710178611e-04


% TW raw data, predhor = 24
0.008320213922550
  -0.031462968682613
  -0.119899769054043
   1.038445160069286
   0.186287938929370

sigma = 0.882584990974591
ACF_TW(2) = -6.403390829217158e-05

%% WHW 2006, prediction horizon 60
% RGS raw data, predhor = 60
  -0.041818188436377
   0.014135156525734
  -0.065718411708037
   0.715325133405660
  -4.189685227801349
   
sigma = 78.767906640296133

ACF_RGS(2) = -5.117865551199228e-04;
bounds = 0.004376923020075;

% TA raw data, predhor = 60
  -0.001494595719600
  -0.070427263542744
  -0.103268450647839
   1.069911163385608
   0.140030683527358

sigma = 0.983962252901425
ACF_TA(2) = 2.340666710178611e-04


% TW raw data, predhor = 60
   0.006455966258724
  -0.043847014073389
  -0.144868533556069
   1.092427188330905
   0.171042158430988
    
sigma = 0.860403161598801
ACF_TW(2) = -6.403390829217158e-05










