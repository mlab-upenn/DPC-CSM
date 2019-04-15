% comp_data.m
% brings prediction of every 12 hours to prediction of every 1 hour form

clc; clear all; close all

weather = 'MSM2006';

%%
load(['TAP_' weather '.mat']);      % load vector with realizations of RGS, 8760 x 1
load(['TWP_' weather '.mat']);  % load raw predictions: 730 x 73, first column useless
load(['TA_' weather '.mat']);      % load vector with realizations of RGS, 8760 x 1
load(['TW_' weather '.mat']);  % load raw predictions: 730 x 73, first column useless
load(['V_' weather '.mat']);      % load vector with realizations of RGS, 8760 x 1

%%
predHor = 60+1;      % used to collect data, must be smaller than (72-12) hours (= 2.5 days)

% bring realizations into "data" form
dataV = nan(length(V)-predHor+1,predHor);
dataTA = nan(length(V)-predHor+1,predHor);
dataTW = nan(length(V)-predHor+1,predHor);
for ii = 1 : length(V)-predHor+1
    dataV(ii,:) = V(ii:ii+predHor-1);
    dataTA(ii,:) = TA(ii:ii+predHor-1);
    dataTW(ii,:) = TW(ii:ii+predHor-1);
end

% save(['dataV_' weather '.mat'],'dataV');
% save(['dataTA_' weather '.mat'],'dataTA');
% save(['dataTW_' weather '.mat'],'dataTW');


%% bring predictions into "data" form
dataTAP = nan(8760,predHor);
dataTWP = nan(8760,predHor);
for ii = 1 : length(V)
    if mod(ii,12) == 0
        dataTAP(ii,:) = TAP(ii/12,12:12+predHor-1);
        dataTWP(ii,:) = TWP(ii/12,12:12+predHor-1);
    else
        dataTAP(ii,:) = TAP(1+floor(ii/12),mod(ii,12):+mod(ii,12)+predHor-1);
        dataTWP(ii,:) = TWP(1+floor(ii/12),mod(ii,12):+mod(ii,12)+predHor-1);
    end
end

% save(['dataTAP_' weather '.mat'],'dataTAP');
% save(['dataTWP_' weather '.mat'],'dataTWP');

%% create dataErr***_KF_***
clc; clear all
weather = 'MSM2006';
load(['dataTAP_KF_' weather '.mat']);
load(['dataTWP_KF_' weather '.mat']);
load(['dataRGSP_KF_' weather '.mat']);
load(['dataTA_' weather '.mat']);
load(['dataTW_' weather '.mat']);
load(['dataV_' weather '.mat']);

dataRGSP_KF(8701:end,:) = []; % vector too long
dataTAP_KF(8701:end,:) = []; % vector too long
dataTWP_KF(8701:end,:) = []; % vector too long

dataErrRGS_KF = dataV - dataRGSP_KF;
dataErrTA_KF = dataTA - dataTAP_KF;
dataErrTW_KF = dataTW - dataTWP_KF;

% save(['dataErrRGS_KF_' weather '.mat'],'dataErrRGS_KF');
% save(['dataErrTA_KF_' weather '.mat'],'dataErrTA_KF');
% save(['dataErrTW_KF_' weather '.mat'],'dataErrTW_KF');


%% create dataErr_***
clc; clear all
weather = 'MSM2006';
load(['dataTAP_' weather '.mat']);
load(['dataTWP_' weather '.mat']);
load(['dataRGSP_' weather '.mat']);
load(['dataTA_' weather '.mat']);
load(['dataTW_' weather '.mat']);
load(['dataV_' weather '.mat']);

dataRGSP(8701:end,:) = []; % vector too long
dataTAP(8701:end,:) = []; % vector too long
dataTWP(8701:end,:) = []; % vector too long

dataErrRGS = dataV - dataRGSP;
dataErrTA = dataTA - dataTAP;
dataErrTW = dataTW - dataTWP;

save(['dataErrRGS_' weather '.mat'],'dataErrRGS');
save(['dataErrTA_' weather '.mat'],'dataErrTA');
save(['dataErrTW_' weather '.mat'],'dataErrTW');