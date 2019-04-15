% fit_AR.m
% file computes fitted AR model for KF, AR2, ...
% used only to fit orignal data, KF-data, AR2-data, not for finding AR2

clc; clear all; close all

weather = 'WHW2006';

%% fitting for SMPC augmentation
load(['dataErrRGS_' weather '.mat']);
load(['dataErrTA_' weather '.mat']);
load(['dataErrTW_' weather '.mat']);
% load(['dataErrRGS_KF_' weather '.mat']);
% load(['dataErrTA_KF_' weather '.mat']);
% load(['dataErrTW_KF_' weather '.mat']);

predHor = 60;
errorRGS = []; 
errorTA = []; 
errorTW = [];
% errorRGS_KF = []; 
% errorTA_KF = []; 
% errorTW_KF = [];
for ii = 1 : size(dataErrTA,1)     % collect "filtered errors"
%     errorRGS_KF = [errorRGS_KF , dataErrRGS_KF(ii,1:predHor)];
%     errorTA_KF = [errorTA_KF , dataErrTA_KF(ii,1:predHor)];
%     errorTW_KF = [errorTW_KF , dataErrTW_KF(ii,1:predHor)];
    errorRGS = [errorRGS , dataErrRGS(ii,1:predHor)];
    errorTA = [errorTA , dataErrTA(ii,1:predHor)];
    errorTW = [errorTW , dataErrTW(ii,1:predHor)];
end

%% AR1_model:
numPrev = 4;

errorRGS_cur = errorRGS(numPrev+1:end)';
X_RGS = [];
for i = 1 : numPrev     % [prev3, prev2, prev1]
    X_RGS = [X_RGS errorRGS(i:end-numPrev+i-1)'];
end

errorTA_cur = errorTA(numPrev+1:end)';
X_TA = [];
for i = 1 : numPrev     % [prev3, prev2, prev1]
    X_TA = [X_TA errorTA(i:end-numPrev+i-1)'];
end

errorTW_cur = errorTW(numPrev+1:end)';
X_TW = [];
for i = 1 : numPrev     % [prev3, prev2, prev1]
    X_TW = [X_TW errorTW(i:end-numPrev+i-1)'];
end

% errorRGS_KF_cur = errorRGS_KF(numPrev+1:end)';
% X_RGS_KF = [];
% for i = 1 : numPrev     % [prev3, prev2, prev1]
%     X_RGS_KF = [X_RGS_KF errorRGS_KF(i:end-numPrev+i-1)'];
% end
% 
% errorTA_KF_cur = errorTA_KF(numPrev+1:end)';
% X_TA_KF = [];
% for i = 1 : numPrev     % [prev3, prev2, prev1]
%     X_TA_KF = [X_TA_KF errorTA_KF(i:end-numPrev+i-1)'];
% end
% 
% errorTW_KF_cur = errorTW_KF(numPrev+1:end)';
% X_TW_KF = [];
% for i = 1 : numPrev     % [prev3, prev2, prev1]
%     X_TW_KF = [X_TW_KF errorTW_KF(i:end-numPrev+i-1)'];
% end

const = ones(length(errorRGS_cur),1);

%%
[b_RGS, bint_RGS, r_RGS, rint_RGS, stats_RGS] = regress(errorRGS_cur,[X_RGS const]);
[b_TA, bint_TA, r_TA, rint_TA, stats_TA] = regress(errorTA_cur,[X_TA const]);
[b_TW, bint_TW, r_TW, rint_TW, stats_TW] = regress(errorTW_cur,[X_TW const]);
[ACF_RGS, lags_RGS, bounds_RGS] = autocorr(r_RGS,length(r_RGS)-1);
[ACF_TA, lags_TA, bounds_TA] = autocorr(r_TA,length(r_TA)-1);
[ACF_TW, lags_TW, bounds_TW] = autocorr(r_TW,length(r_TW)-1);

% [b_RGS_KF, bint_RGS_KF, r_RGS_KF, rint_RGS_KF, stats_RGS_KF] = regress(errorRGS_KF_cur,[X_RGS_KF const]);
% [b_TA_KF, bint_TA_KF, r_TA_KF, rint_TA_KF, stats_TA_KF] = regress(errorTA_KF_cur,[X_TA_KF const]);
% [b_TW_KF, bint_TW_KF, r_TW_KF, rint_TW_KF, stats_TW_KF] = regress(errorTW_KF_cur,[X_TW_KF const]);
% [ACF_RGS_KF, lags_RGS_KF, bounds_RGS_KF] = autocorr(r_RGS_KF,length(r_RGS_KF)-1);
% [ACF_TA_KF, lags_TA_KF, bounds_TA_KF] = autocorr(r_TA_KF,length(r_TA_KF)-1);
% [ACF_TW_KF, lags_TW_KF, bounds_TW_KF] = autocorr(r_TW_KF,length(r_TW_KF)-1);

% % save residua for samples
% errDeltaRGS = r_RGS;
% save(['errDeltaRGS_' num2str(predHor) '_' weather],'errDeltaRGS');
% errDeltaTA = r_TA;
% save(['errDeltaTA_' num2str(predHor) '_' weather],'errDeltaTA');
% errDeltaTW = r_TW;
% save(['errDeltaTW_' num2str(predHor) '_' weather],'errDeltaTW');