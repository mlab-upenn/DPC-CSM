%% noisyMeasurement.m
% Date: 2. August 2012
% Author: Xiaojing Zhang
% creates noisy measurement data based on "real" data retrieved from OCWDB
% only interested in 2007 data since we don't need it for 2006
% noisy measurement fuer Kalman Filter (KF) estimation of error

clear all
close all
clc

weather = 'WHW2007';

%% 
cd(weather);
load(['TAreal_' weather '.mat'])
load(['TWreal_' weather '.mat'])
load(['Vreal_' weather '.mat'])

TAmeas = TA + 0.1*randn(length(TA),1);      % 2*sigma = 0.2
TWmeas = TW + 0.1*randn(length(TW),1);
Vmeas = V + 40*randn(length(V),1);          % sigma = 40 W/m

Vmeas = max(Vmeas,0);   % never smaller than 0
for ii = 1 : length(V)
    if V(ii) < 0.5      % night
        Vmeas(ii) = 0;
    end
end


figure(); hold on;
plot(TA,'b','LineWidth',2)
plot(TAmeas,':r','LineWidth',2)
legend('real TA','measured TA')

figure(); hold on;
plot(TW,'b','LineWidth',2)
plot(TWmeas,':r','LineWidth',2)
legend('real TW','measured TW')

figure(); hold on;
plot(V,'b','LineWidth',2)
plot(Vmeas,':r','LineWidth',2)
legend('real RG','measured RG')

save(['TAmeas_' weather '.mat'],'TAmeas')
save(['TWmeas_' weather '.mat'],'TWmeas')
save(['Vmeas_' weather '.mat'],'Vmeas')

cd ..