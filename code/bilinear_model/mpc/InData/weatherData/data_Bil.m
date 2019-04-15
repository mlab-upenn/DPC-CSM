%% data_Bil.m
% constructs dataRGSP_Bil, dataTAP_Bil, dataTWP_Bil for use in bilinearity of system
% takes predicted data, fits error model and improves predicted data
% - uses measurements

% Author: Xiaojing George Zhang
% Date: 13. Aug 2012


clear all; 
close all; 
clc;

weather  = 'WHW2007';
predhor = 60;

%% load error, predicted values and 
cd(weather)
load(['dataRGSP_' weather '.mat'])
load(['dataTAP_' weather '.mat'])
load(['dataTWP_' weather '.mat'])

load(['dataRGSP_KF_' weather '.mat'])
load(['dataTAP_KF_' weather '.mat'])
load(['dataTWP_KF_' weather '.mat'])

load(['dataErrTA_' weather '.mat'])
load(['dataErrTW_' weather '.mat'])
load(['dataErrRGS_' weather '.mat'])

load(['errRGS_measKF_' weather '.mat'])
load(['errTA_measKF_' weather '.mat'])
load(['errTW_measKF_' weather '.mat'])

load(['dataVreal_' weather '.mat'])
load(['dataTAreal_' weather '.mat'])
load(['dataTWreal_' weather '.mat'])

cd ..
% cut off predictions to 8700x61
dataRGSP(8701:end,:) = [];
dataTAP(8701:end,:) = [];
dataTWP(8701:end,:) = [];

% dataRGSP_KF(8701:end,:) = [];
% dataTAP_KF(8701:end,:) = [];
% dataTWP_KF(8701:end,:) = [];

%% build dataTAP_Bil, dataTWP_Bil, dataRGSP_Bil
% TA MSM 2007: AR model of order 4
if strcmp(weather,'MSM2007')
    if predhor == 24
        b4_TA = -0.005798218875725;
        b3_TA = -0.073334909041129;
        b2_TA = -0.185413392371912;
        b1_TA = 1.116897869277878;
        b0_TA = 0.026546339096447;  
    elseif predhor == 60
        b4_TA = 0.005036885519147;
        b3_TA = -0.074657644688327;
        b2_TA = -0.203705866295020;
        b1_TA = 1.143846996362533;
        b0_TA = 0.030625458333543; 
    else
        error('do not recognize predhor')
    end
elseif strcmp(weather,'WHW2007')
    if predhor == 24
        b4_TA = -0.006658251730190;
        b3_TA = -0.062359836733178;
        b2_TA = -0.088413137437204;
        b1_TA = 1.037095643071405;
        b0_TA = 0.168355788215680;
    elseif predhor == 60
        b4_TA = -0.001494595719600;
        b3_TA = -0.070427263542744;
        b2_TA = -0.103268450647839;
        b1_TA = 1.069911163385608;
        b0_TA = 0.168355788215680;
    else
        error('do not recognize predhor')
    end
else
    error('Do not recognize weather')
end
% treat dataTAP_Bil
dataTAP_KF_AR4 = zeros(size(dataTA));
err_prev1 = 0;       % uses previous error to get better predictions
err_prev2 = 0;
err_prev3 = 0;
err_prev4 = 0;
for ii = 1 : size(dataTAP_KF_AR4,1)
    for jj = 1 : size(dataTAP_KF_AR4,2)
%         err_prev = b1_TA*err_prev + b0_TA;      % update
        err_update = b1_TA*err_prev1 + b2_TA*err_prev2 + b3_TA*err_prev3 + b4_TA*err_prev4 + b0_TA;
        dataTAP_KF_AR4(ii,jj) = dataTAP(ii,jj) + err_update;
        err_prev4 = err_prev3;
        err_prev3 = err_prev2;
        err_prev2 = err_prev1;
        err_prev1 = err_update;
    end
    err_prev4 = err_prev3;
    err_prev3 = err_prev2;
    err_prev2 = err_prev1;
%     err_prev1 = dataErrTA(ii,1);     % update error for next hour
    err_prev1 = errTA_measKF(ii);
end

figure(); hold on;
title('TA first step')
plot(dataTA(:,1),'b','LineWidth',2)
plot(dataTAP(:,1),'--g','LineWidth',2)
plot(dataTAP_KF_AR4(:,1),'-.r','LineWidth',2)
plot(dataTAP_KF(:,1),':c','LineWidth',2)
% legend('TA Realization','TA Prediction','TA Prediction Bil','TA Prediction KF')
legend('TA Realization','raw TA Prediction','TA Prediction KF AR4','TA Prediction Vanessa KF')
hold off

jj = 15;
figure(); hold on;
title('TA first step')
plot(dataTA(jj,:),'b','LineWidth',2)
plot(dataTAP(jj,:),'--g','LineWidth',2)
plot(dataTAP_KF_AR4(jj,:),'-.r','LineWidth',2)
plot(dataTAP_KF(jj,:),':c','LineWidth',2)
% legend('TA Realization','TA Prediction','TA Prediction Bil','TA Prediction KF')
legend('TA Realization','raw TA Prediction','TA Prediction KF AR4','TA Prediction KF')
hold off


%% build dataTAP_Bil, dataTWP_Bil, dataRGSP_Bil
% TA MSM 2007: AR model of order 4
if strcmp(weather, 'MSM2007')
    if predhor == 24
        b4_TW = -0.021789228921885;
        b3_TW = -0.026285133043742;
        b2_TW = -0.177805832504286;
        b1_TW = 1.127068855353699;
        b0_TW = 0.339531486192709; 
    elseif predhor == 60
        b4_TW = -0.016161958287229;
        b3_TW = -0.032380350240663;
        b2_TW = -0.195734305981573;
        b1_TW = 1.155754858223748;
        b0_TW = 0.314040430720039; 
    else
        error('do not recognize predhor')
    end
elseif strcmp(weather, 'WHW2007')
    if predhor == 24
        b4_TW = 0.008320213922550;
        b3_TW = -0.031462968682613;
        b2_TW = -0.119899769054043;
        b1_TW = 1.038445160069286;
        b0_TW = 0.186287938929370;
    elseif predhor == 60
        b4_TW = 0.006455966258724;
        b3_TW = -0.043847014073389;
        b2_TW = -0.144868533556069;
        b1_TW = 1.092427188330905;
        b0_TW = 0.171042158430988;
    else
        error('do not recognize predhor')
    end
else
    error('do not recognize weather')
end
 

% treat dataTAP_Bil
dataTWP_KF_AR4 = zeros(size(dataTW));
err_prev1 = 0;       % uses previous error to get better predictions
err_prev2 = 0;
err_prev3 = 0;
err_prev4 = 0;
for ii = 1 : size(dataTWP_KF_AR4,1)
    for jj = 1 : size(dataTWP_KF_AR4,2)
%         err_prev = b1_TA*err_prev + b0_TA;      % update
        err_update = b1_TW*err_prev1 + b2_TW*err_prev2 + b3_TW*err_prev3 + b4_TW*err_prev4 + b0_TW;
        dataTWP_KF_AR4(ii,jj) = dataTWP(ii,jj) + err_update;
        err_prev4 = err_prev3;
        err_prev3 = err_prev2;
        err_prev2 = err_prev1;
        err_prev1 = err_update;
    end
    err_prev4 = err_prev3;
    err_prev3 = err_prev2;
    err_prev2 = err_prev1;
    err_prev1 = errTW_measKF(ii,1);     % update error for next hour
end

figure(); hold on;
title('TW first step')
plot(dataTW(:,1),'b','LineWidth',2)
plot(dataTWP(:,1),'--g','LineWidth',2)
plot(dataTWP_KF_AR4(:,1),'-.r','LineWidth',2)
plot(dataTWP_KF(:,1),':c','LineWidth',2)
% legend('TW Realization','TW Prediction','TW Prediction Bil','TW Prediction KF')
legend('TW Realization','raw TW prediction','TW Prediction KF AR4', 'TW Predictin Vanessa KF')
hold off

jj = 9;
figure(); hold on;
title('TW first step')
plot(dataTW(jj,:),'b','LineWidth',2)
plot(dataTWP(jj,:),'--g','LineWidth',2)
plot(dataTWP_KF_AR4(jj,:),'-.r','LineWidth',2)
plot(dataTWP_KF(jj,:),':c','LineWidth',2)
legend('TW Realization','raw TW Prediction','TA Prediction KF AR4')
hold off


%% build dataTAP_Bil, dataTWP_Bil, dataRGSP_Bil
% TA MSM 2007: AR model of order 4
if strcmp(weather, 'MSM2007')
    if predhor == 24
        b4_RGS = -0.085522186754388;
        b3_RGS = -0.006576114457115;
        b2_RGS = -0.110220643581048;
        b1_RGS = 0.791052360508261;
        b0_RGS = 2.595787479120149;  
    elseif predhor == 60
        b4_RGS = -0.092933939337812;
        b3_RGS = 0.012435958044596;
        b2_RGS = -0.132082680832671;
        b1_RGS = 0.834136286496958;
        b0_RGS = 2.131798066846339;
    else
        error('do not recognize predhor')
    end
elseif strcmp(weather, 'WHW2007')
    if predhor == 24
        b4_RGS = -0.045264917827551;
        b3_RGS = 0.023166328048953;
        b2_RGS = -0.068075467934311;
        b1_RGS = 0.685448697438124;
        b0_RGS = -4.700563969468201;
    elseif predhor == 60
        b4_RGS = -0.041818188436377;
        b3_RGS = 0.014135156525734;
        b2_RGS = -0.065718411708037;
        b1_RGS = 0.715325133405660;
        b0_RGS = -4.189685227801349;
    else
        error('do not recognize predhor')
    end
else
    error('do not recognize weather')
end


% treat dataTAP_Bil
dataRGSP_KF_AR4 = zeros(size(dataV));
err_prev1 = 0;       % uses previous error to get better predictions
err_prev2 = 0;
err_prev3 = 0;
err_prev4 = 0;
for ii = 1 : size(dataRGSP_KF_AR4,1)
    for jj = 1 : size(dataRGSP_KF_AR4,2)
%         err_prev = b1_TA*err_prev + b0_TA;      % update
        err_update = b1_RGS*err_prev1 + b2_RGS*err_prev2 + b3_RGS*err_prev3 + b4_RGS*err_prev4 + b0_RGS;
        dataRGSP_KF_AR4(ii,jj) = dataRGSP(ii,jj) + err_update;
        err_prev4 = err_prev3;
        err_prev3 = err_prev2;
        err_prev2 = err_prev1;
        err_prev1 = err_update;
    end
    err_prev4 = err_prev3;
    err_prev3 = err_prev2;
    err_prev2 = err_prev1;
    err_prev1 = errRGS_measKF(ii,1);     % update error for next hour
end

%RGS only positive
dataRGSP_KF_AR4 = max(dataRGSP_KF_AR4,0);

%RGSP = 0 for night
for ii = 1 : size(dataV,1)
    for jj = 1 : size(dataV,2)
        if abs(dataV(ii,jj)) <= 0.5
            dataRGSP_KF_AR4(ii,jj) = 0;
        end
    end
end

figure(); hold on;
title('RGS first step')
plot(dataV(:,1),'b','LineWidth',2)
plot(dataRGSP(:,1),'--g','LineWidth',2)
plot(dataRGSP_KF_AR4(:,1),'-.r','LineWidth',2)
plot(dataRGSP_KF(:,1),':c','LineWidth',2)
legend('RGS Realization','RGS raw Prediction','RGS Prediction KF AR4','RGS Prediction Vanessa KF')
hold off

jj = 9;
figure(); hold on;
title('RGS first step')
plot(dataV(jj,:),'b','LineWidth',2)
plot(dataRGSP(jj,:),'--g','LineWidth',2)
plot(dataRGSP_KF_AR4(jj,:),'-.r','LineWidth',2)
plot(dataRGSP_KF(jj,:),':c','LineWidth',2)
legend('RGS Realization','RGS Prediction','RGS Prediction KF AR4')
hold off

%%
cd(weather)
% save(['dataRGSP_KF_AR4_' num2str(predhor) '_' weather '.mat'], 'dataRGSP_KF_AR4')
% save(['dataTAP_KF_AR4_' num2str(predhor) '_' weather '.mat'], 'dataTAP_KF_AR4')
% save(['dataTWP_KF_AR4_' num2str(predhor) '_' weather '.mat'], 'dataTWP_KF_AR4')
cd ..


%% make plots for Report
jj = 15
figure(); hold on;
title('Air Temperature (TA): Realization and Forecast', 'FontSize',18)
plot(dataTA(jj,:),'b','LineWidth',2)
plot(dataTAP(jj,:),'-.k','LineWidth',2)
h_legend = legend('Realization $v_{t+i}$','Forecast $\bar{v}_{t+i|t}$');
set(h_legend,'FontSize',16,'Interpreter','LATEX')
xlabel('time steps [h]','FontSize',18)
ylabel('temperature [degC]','FontSize',18)
hold off

figure(); hold on;
title('Air Temperature (TA): Realization, Forecast and Prediction', 'FontSize',18)
plot(dataTA(jj,:),'b','LineWidth',2)
plot(dataTAP(jj,:),'-.k','LineWidth',2)
plot(dataTAP_KF_AR4(jj,:),'--r','LineWidth',2)
h_legend = legend('Realization $v_{t+i}$','Forecast $\bar{v}_{t+i|t}$',...
                    'Prediction $\hat{v}_{t+i|t}$');
set(h_legend,'FontSize',16,'Interpreter','LATEX')
xlabel('time steps [h]','FontSize',18)
ylabel('temperature [degC]','FontSize',18)
hold off

figure(); hold on;
title('Air Temperature (TA): One-step Prediction', 'FontSize',18)
plot(dataTA(jj,:),'b','LineWidth',2)
plot(dataTAP_KF_AR4(jj,:),'--r','LineWidth',2)
plot(dataTAP_KF_AR4(jj:jj+60,1),'-.m','LineWidth',2)
h_legend = legend('Realization $v_{t+i}$', 'Prediction $\hat{v}_{t+i|t}$', ...
                    'Prediction $\hat{v}_{t+i+1|t+i}$');
set(h_legend,'FontSize',16,'Interpreter','LATEX')
xlabel('time steps [h]','FontSize',18)
ylabel('temperature [degC]','FontSize',18)
hold off

