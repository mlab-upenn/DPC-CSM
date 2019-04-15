%% data_Bil.m
% constructs dataRGSP_Bil, dataTAP_Bil, dataTWP_Bil for use in bilinearity of system
% takes predicted data, fits error model and improves predicted data

% Author: Xiaojing George Zhang
% Date: 31. May 2012


clear all; close all; clc;

predhor = 24;
%% load error, predicted values and 

load('dataRGSP_MSM2007.mat')
load('dataTAP_MSM2007.mat')
load('dataTWP_MSM2007.mat')

load('dataRGSP_KF_MSM2007.mat')
load('dataTAP_KF_MSM2007.mat')
load('dataTWP_KF_MSM2007.mat')

load('dataErrTA_MSM2007.mat')
load('dataErrTW_MSM2007.mat')
load('dataErrRGS_MSM2007.mat')

load('dataV_MSM2007.mat')
load('dataTA_MSM2007.mat')
load('dataTW_MSM2007.mat')


% cut off predictions to 8700x61
dataRGSP(8701:end,:) = [];
dataTAP(8701:end,:) = [];
dataTWP(8701:end,:) = [];

dataRGSP_KF(8701:end,:) = [];
dataTAP_KF(8701:end,:) = [];
dataTWP_KF(8701:end,:) = [];

%% build dataTAP_Bil, dataTWP_Bil, dataRGSP_Bil
% TA MSM 2007: AR model of order 4
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
    error('do not recognize prediction horizon')
end
% treat dataTAP_Bil
dataTAP_AR4 = zeros(size(dataTA));
err_prev1 = 0;       % uses previous error to get better predictions
err_prev2 = 0;
err_prev3 = 0;
err_prev4 = 0;
for ii = 1 : size(dataTAP_AR4,1)
    for jj = 1 : size(dataTAP_AR4,2)
%         err_prev = b1_TA*err_prev + b0_TA;      % update
        err_update = b1_TA*err_prev1 + b2_TA*err_prev2 + b3_TA*err_prev3 + b4_TA*err_prev4 + b0_TA;
        dataTAP_AR4(ii,jj) = dataTAP(ii,jj) + err_update;
        err_prev4 = err_prev3;
        err_prev3 = err_prev2;
        err_prev2 = err_prev1;
        err_prev1 = err_update;
    end
    err_prev4 = err_prev3;
    err_prev3 = err_prev2;
    err_prev2 = err_prev1;
    err_prev1 = dataErrTA(ii,1);     % update error for next hour
end

figure(); hold on;
title('TA first step')
plot(dataTA(:,1),'b','LineWidth',2)
% plot(dataTAP(:,1),'--g','LineWidth',2)
plot(dataTAP_AR4(:,1),'-.r','LineWidth',2)
plot(dataTAP_KF(:,1),':c','LineWidth',2)
% legend('TA Realization','TA Prediction','TA Prediction Bil','TA Prediction KF')
legend('TA Realization','TA Prediction Bil','TA Prediction KF')
hold off

jj = 15;
figure(); hold on;
title('TA first step')
plot(dataTA(jj,:),'b','LineWidth',2)
% plot(dataTAP(jj,:),'--r','LineWidth',2)
plot(dataTAP_AR4(jj,:),'-.r','LineWidth',2)
plot(dataTAP_KF(jj,:),':c','LineWidth',2)
% legend('TA Realization','TA Prediction','TA Prediction Bil','TA Prediction KF')
legend('TA Realization','TA Prediction Bil','TA Prediction KF')
hold off


%% build dataTAP_Bil, dataTWP_Bil, dataRGSP_Bil
% TA MSM 2007: AR model of order 4
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
% treat dataTAP_Bil
dataTWP_AR4 = zeros(size(dataTW));
err_prev1 = 0;       % uses previous error to get better predictions
err_prev2 = 0;
err_prev3 = 0;
err_prev4 = 0;
for ii = 1 : size(dataTWP_AR4,1)
    for jj = 1 : size(dataTWP_AR4,2)
%         err_prev = b1_TA*err_prev + b0_TA;      % update
        err_update = b1_TW*err_prev1 + b2_TW*err_prev2 + b3_TW*err_prev3 + b4_TW*err_prev4 + b0_TW;
        dataTWP_AR4(ii,jj) = dataTWP(ii,jj) + err_update;
        err_prev4 = err_prev3;
        err_prev3 = err_prev2;
        err_prev2 = err_prev1;
        err_prev1 = err_update;
    end
    err_prev4 = err_prev3;
    err_prev3 = err_prev2;
    err_prev2 = err_prev1;
    err_prev1 = dataErrTW(ii,1);     % update error for next hour
end

figure(); hold on;
title('TW first step')
plot(dataTW(:,1),'b','LineWidth',2)
% plot(dataTWP(:,1),'--r','LineWidth',2)
plot(dataTWP_AR4(:,1),'-.r','LineWidth',2)
plot(dataTWP_KF(:,1),':c','LineWidth',2)
% legend('TW Realization','TW Prediction','TW Prediction Bil','TW Prediction KF')
legend('TW Realization','TW Prediction Bil','TW Prediction KF')
hold off

jj = 9;
figure(); hold on;
title('TW first step')
plot(dataTW(jj,:),'b')
plot(dataTWP(jj,:),'--r')
plot(dataTWP_AR4(jj,:),'-.g')
legend('TW Realization','TW Prediction','TA Prediction Bil')
hold off

%% build dataTAP_Bil, dataTWP_Bil, dataRGSP_Bil
% TA MSM 2007: AR model of order 4
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
    error('do not recognize pred hor')
end
% treat dataTAP_Bil
dataRGSP_AR4 = zeros(size(dataV));
err_prev1 = 0;       % uses previous error to get better predictions
err_prev2 = 0;
err_prev3 = 0;
err_prev4 = 0;
for ii = 1 : size(dataRGSP_AR4,1)
    for jj = 1 : size(dataRGSP_AR4,2)
%         err_prev = b1_TA*err_prev + b0_TA;      % update
        err_update = b1_RGS*err_prev1 + b2_RGS*err_prev2 + b3_RGS*err_prev3 + b4_RGS*err_prev4 + b0_RGS;
        dataRGSP_AR4(ii,jj) = dataRGSP(ii,jj) + err_update;
        err_prev4 = err_prev3;
        err_prev3 = err_prev2;
        err_prev2 = err_prev1;
        err_prev1 = err_update;
    end
    err_prev4 = err_prev3;
    err_prev3 = err_prev2;
    err_prev2 = err_prev1;
    err_prev1 = dataErrRGS(ii,1);     % update error for next hour
end

%RGS only positive
dataRGSP_AR4 = max(dataRGSP_AR4,0);

%RGSP = 0 for night
for ii = 1 : size(dataV,1)
    for jj = 1 : size(dataV,2)
        if abs(dataV(ii,jj)) <= 0.5
            dataRGSP_AR4(ii,jj) = 0;
        end
    end
end

figure(); hold on;
title('RGS first step')
plot(dataV(:,1),'b','LineWidth',2)
% plot(dataRGSP(:,1),'--r','LineWidth',2)
plot(dataRGSP_AR4(:,1),'-.r','LineWidth',2)
plot(dataRGSP_KF(:,1),':c','LineWidth',2)
% legend('RGS Realization','RGS Prediction','RGS Prediction Bil','RGS Prediction KF')
legend('RGS Realization','RGS Prediction Bil','RGS Prediction KF')
hold off

jj = 9;
figure(); hold on;
title('RGS first step')
plot(dataV(jj,:),'b')
plot(dataRGSP(jj,:),'--r')
plot(dataRGSP_AR4(jj,:),'-.g')
legend('RGS Realization','RGS Prediction','RGS Prediction Bil')
hold off

%%
save(['dataRGSP_AR4_' num2str(predhor) '_MSM2007.mat'], 'dataRGSP_AR4')
save(['dataTAP_AR4_' num2str(predhor) '_MSM2007.mat'], 'dataTAP_AR4')
save(['dataTWP_AR4_' num2str(predhor) '_MSM2007.mat'], 'dataTWP_AR4')


%% compute error RGSP
% start midnight
err_Bil = [];
err_orig = [];
for ii = 1 : 24 : size(dataTWP_AR4,1)       % starting at night
    err_Bil = [err_Bil
           dataV(ii,:) - dataRGSP_AR4(ii,:)   ];
    err_orig = [err_orig
           dataErrRGS(ii,:)   ];
end
err_Bil_mean = mean(err_Bil);
err_orig_mean = mean(err_orig);
mean(err_Bil_mean)
mean(err_orig_mean)

figure(); hold on;
title('prediction error RGS from midnight')
for ii = 1 : size(err_Bil,1)
   plot(err_orig(ii,:),'b')
   plot(err_Bil(ii,:),'r')
end
hold off

figure(); hold on;
title('mean prediction error RGS from midnight')
plot(err_orig_mean,'b')
plot(err_Bil_mean,'r')
legend('unfiltered prediction error', 'AR1-filtered prediction error')
hold off

% compute error
err_Bil = [];
err_orig = [];
for ii = 13 : 24 : size(dataTWP_AR4,1)      % starting noon
    err_Bil = [err_Bil
           dataV(ii,:) - dataRGSP_AR4(ii,:)   ];
    err_orig = [err_orig
           dataErrRGS(ii,:)   ];
end
err_Bil_mean = mean(err_Bil);
err_orig_mean = mean(err_orig);
mean(err_Bil_mean)
mean(err_orig_mean)

figure(); hold on;
title('prediction error RGS from midnight')
for ii = 1 : size(err_Bil,1)
   plot(err_orig(ii,:),'b')
   plot(err_Bil(ii,:),'r')
end
hold off

figure(); hold on;
title('mean prediction error RGS from noon')
plot(err_orig_mean,'b')
plot(err_Bil_mean,'r')
legend('unfiltered prediction error', 'AR1-filtered prediction error')
hold off

%% compute error TA
% start midnight
err_Bil = [];
err_orig = [];
for ii = 1 : 24 : size(dataTAP_AR4,1)       % starting at night
    err_Bil = [err_Bil
           dataTA(ii,:) - dataTAP_AR4(ii,:)   ];
    err_orig = [err_orig
           dataErrTA(ii,:)   ];
end
err_Bil_mean = mean(err_Bil);
err_orig_mean = mean(err_orig);
mean(err_Bil_mean)
mean(err_orig_mean)

figure(); hold on;
title('prediction error TA from midnight')
for ii = 1 : size(err_Bil,1)
   plot(err_orig(ii,:),'b')
   plot(err_Bil(ii,:),'r')
end
hold off

figure(); hold on;
title('mean prediction error TA from midnight')
plot(err_orig_mean,'b')
plot(err_Bil_mean,'r')
legend('unfiltered prediction error', 'AR1-filtered prediction error')
hold off

% compute error
err_Bil = [];
err_orig = [];
for ii = 13 : 24 : size(dataTAP_AR4,1)      % starting noon
    err_Bil = [err_Bil
           dataTA(ii,:) - dataTAP_AR4(ii,:)   ];
    err_orig = [err_orig
           dataErrTA(ii,:)   ];
end
err_Bil_mean = mean(err_Bil);
err_orig_mean = mean(err_orig);
mean(err_Bil_mean)
mean(err_orig_mean)

figure(); hold on;
title('prediction error TA from midnight')
for ii = 1 : size(err_Bil,1)
   plot(err_orig(ii,:),'b')
   plot(err_Bil(ii,:),'r')
end
hold off

figure(); hold on;
title('mean prediction error TA from noon')
plot(err_orig_mean,'b')
plot(err_Bil_mean,'r')
legend('unfiltered prediction error', 'AR1-filtered prediction error')
hold off

%% compute error TW
% start midnight
err_Bil = [];
err_orig = [];
for ii = 1 : 24 : size(dataTAP_AR4,1)       % starting at night
    err_Bil = [err_Bil
           dataTW(ii,:) - dataTWP_AR4(ii,:)   ];
    err_orig = [err_orig
           dataErrTW(ii,:)   ];
end
err_Bil_mean = mean(err_Bil);
err_orig_mean = mean(err_orig);
mean(err_Bil_mean)
mean(err_orig_mean)

figure(); hold on;
title('prediction error TW from midnight')
for ii = 1 : size(err_Bil,1)
   plot(err_orig(ii,:),'b')
   plot(err_Bil(ii,:),'r')
end
hold off

figure(); hold on;
title('mean prediction error TW from midnight')
plot(err_orig_mean,'b')
plot(err_Bil_mean,'r')
legend('unfiltered prediction error', 'AR1-filtered prediction error')
hold off

% compute error
err_Bil = [];
err_orig = [];
for ii = 13 : 24 : size(dataTAP_AR4,1)      % starting noon
    err_Bil = [err_Bil
           dataTW(ii,:) - dataTWP_AR4(ii,:)   ];
    err_orig = [err_orig
           dataErrTW(ii,:)   ];
end
err_Bil_mean = mean(err_Bil);
err_orig_mean = mean(err_orig);
mean(err_Bil_mean)
mean(err_orig_mean)

figure(); hold on;
title('prediction error TW from midnight')
for ii = 1 : size(err_Bil,1)
   plot(err_orig(ii,:),'b')
   plot(err_Bil(ii,:),'r')
end
hold off

figure(); hold on;
title('mean prediction error TW from noon')
plot(err_orig_mean,'b')
plot(err_Bil_mean,'r')
legend('unfiltered prediction error', 'AR1-filtered prediction error')
hold off

%%
% save('dataRGSP_Bil_MSM2007.mat', 'dataRGSP_Bil')
% save('dataTAP_Bil_MSM2007.mat', 'dataTAP_Bil')
% save('dataTWP_Bil_MSM2007.mat', 'dataTWP_Bil')


