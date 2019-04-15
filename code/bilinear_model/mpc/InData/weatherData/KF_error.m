%% KF_error.m
% Author: Xiaojng Zhang
% Date: 02. August 2012

clc
close all
clear all

weather = 'MSM2007';

%% load data
cd(weather);
load(['Vreal_' weather '.mat'])
load(['TAreal_' weather '.mat'])
load(['TWreal_' weather '.mat'])

load(['Vmeas_' weather '.mat'])
load(['TAmeas_' weather '.mat'])
load(['TWmeas_' weather '.mat'])

load(['dataRGSP_' weather '.mat'])
load(['dataTAP_' weather '.mat'])
load(['dataTWP_' weather '.mat'])

Vpred = dataRGSP(:,1);
TApred = dataTAP(:,1);
TWpred = dataTWP(:,1);

% figure(); hold on;
% plot(V,'b','LineWidth',2)
% plot(Vmeas,':r','LineWidth',2)
% plot(Vpred,'--g','LineWidth',2)
% legend('V real','V meas','V pred')
% hold off
% 
% figure(); hold on;
% plot(TA,'b','LineWidth',2)
% plot(TAmeas,':r','LineWidth',2)
% plot(TApred,'--g','LineWidth',2)
% legend('TA real','TA meas','TA pred')
% hold off
% 
% figure(); hold on;
% plot(TW,'b','LineWidth',2)
% plot(TWmeas,':r','LineWidth',2)
% plot(TWpred,'--g','LineWidth',2)
% legend('TW real','TW meas','TW pred')
% hold off

%% compute error v_tilde
errRGS_real = V - Vpred;
errRGS_meas = Vmeas - Vpred;

errTA_real = TA - TApred;
errTA_meas = TAmeas - TApred;

errTW_real = TW - TWpred;
errTW_meas = TWmeas - TWpred;

% figure(); hold on;
% plot(errRGS_real,'r','LineWidth',2)
% plot(errRGS_meas,'--g','LineWidth',2)
% legend('Error V real','Error V meas')
% hold off
% 
% figure(); hold on;
% plot(errTA_real,'r','LineWidth',2)
% plot(errTA_meas,'--g','LineWidth',2)
% legend('Error TA real','Error TA meas')
% hold off
% 
% 
% figure(); hold on;
% plot(errTW_real,'r','LineWidth',2)
% plot(errTW_meas,'--g','LineWidth',2)
% legend('Error TW real','Error TW meas')
% hold off

%% 

if strcmp(weather,'MSM2007')
    b4_RGS = -0.085522186754388;
    b3_RGS = -0.006576114457115;
    b2_RGS = -0.110220643581048;
    b1_RGS = 0.791052360508261;
    b0_RGS = 2.595787479120149;  
    std_RGS = 76.498961265170635;   % process noise

    b4_TA = -0.005798218875725;
    b3_TA = -0.073334909041129;
    b2_TA = -0.185413392371912;
    b1_TA = 1.116897869277878;
    b0_TA = 0.026546339096447;  
    std_TA = 0.962565568254634;     % process noise

    b4_TW = -0.021789228921885;
    b3_TW = -0.026285133043742;
    b2_TW = -0.177805832504286;
    b1_TW = 1.127068855353699;
    b0_TW = 0.339531486192709; 
    std_TW = 1.060506242002127;     % process noise
    
elseif strcmp(weather,'WHW2007')
    b4_RGS = -0.045264917827551;
    b3_RGS = 0.023166328048953;
    b2_RGS = -0.068075467934311;
    b1_RGS = 0.685448697438124;
    b0_RGS = -4.700563969468201;
    std_RGS = 77.909579020994158;   % process noise
            
    b4_TA = -0.006658251730190;
    b3_TA = -0.062359836733178;
    b2_TA = -0.088413137437204;
    b1_TA = 1.037095643071405;
    b0_TA = 0.168355788215680;
    std_TA = 0.996752980816940;     % process noise
           
    b4_TW = 0.008320213922550;
    b3_TW = -0.031462968682613;
    b2_TW = -0.119899769054043;
    b1_TW = 1.038445160069286;
    b0_TW = 0.186287938929370;
    std_TW = 0.882584990974591;     % process noise
else
    error('Do not recognize Weather type')
end

%% build KF for TA, TW, and 


% state space model with all three cases
% error+ = F*error + w_bar*u + Q*w, where w~N(0,I), u=1  % real error
% error_meas = C*error + H*w + v, where E(vv) = R
F_RGS = [   b1_RGS b2_RGS b3_RGS  b4_RGS
                 1      0      0       0
                 0      1      0       0
                 0      0      1       0];
F_TA = [     b1_TA  b2_TA  b3_TA  b4_TA
                 1      0      0      0 
                 0      1      0      0
                 0      0      1      0];
F_TW = [     b1_TW  b2_TW  b3_TW  b4_TW 
                     1      0      0      0
                     0      1      0      0
                     0      0      1      0];
F = blkdiag(F_RGS, F_TA, F_TW);
    
wbar_RGS = [b0_RGS ; 0 ; 0 ; 0];
wbar_TA = [b0_TA ; 0 ; 0 ; 0];
wbar_TW = [b0_TW ; 0 ; 0 ; 0];
w_bar = [wbar_RGS ; wbar_TA ; wbar_TW];

Q_RGS = [std_RGS ; 0 ; 0 ; 0];
Q_TA = [std_TA ; 0 ; 0 ; 0];
Q_TW = [std_TW ; 0 ; 0 ; 0];
Q = blkdiag(Q_RGS, Q_TA, Q_TW);
    
C_RGS = [1 0 0 0];
C = blkdiag(C_RGS,C_RGS,C_RGS);
D = zeros(3,1);
H = zeros(3,3);
    
R = diag([40^2,0.1^2,0.1^2]);            % R = E[vv^T]

%% construct KF
    
sys = ss(F,[w_bar Q],C,[],-1);
[kalmf,L,P] = kalman(sys,eye(3),R);
kalmf = kalmf(1:3,:);   % only want first couple of steps
u = ones(1,length(V));
t = 0: 1 : length(V)-1;
x0 = [  errRGS_real(1); 0 ; 0 ; 0
        errTA_real(1); 0 ; 0 ; 0
        errTW_real(1); 0 ; 0 ; 0];
err_est = lsim(kalmf, [u; errRGS_meas' ; errTA_meas' ; errTW_meas'], t, zeros(12,1));
    
%% post process RGS predictions
for ii = 1 : length(V)
    if V(ii) < 0.5
       err_est(ii,1) = 0; 
    end
end
    
figure(); hold on;
plot(errRGS_real,'r','LineWidth',2)
plot(errRGS_meas,'--g','LineWidth',2)
plot(err_est(:,1),':b','LineWidth',2)
legend('Error V real','Error V meas','Error V est')
hold off
    
figure(); hold on;
plot(errTA_real,'r','LineWidth',2)
plot(errTA_meas,'--g','LineWidth',2)
plot(err_est(:,2),':b','LineWidth',2)
legend('Error TA real','Error TA meas','Error V est')
hold off

figure(); hold on;
plot(errTW_real,'r','LineWidth',2)
plot(errTW_meas,'--g','LineWidth',2)
plot(err_est(:,3),':b','LineWidth',2)
legend('Error TW real','Error TW meas','Error V est')
hold off
    
    
errRGS_measKF = err_est(:,1);
errTA_measKF = err_est(:,2);
errTW_measKF = err_est(:,3);

% save(['errRGS_measKF_' weather '.mat'],'errRGS_measKF');
% save(['errTA_measKF_' weather '.mat'],'errTA_measKF');
% save(['errTW_measKF_' weather '.mat'],'errTW_measKF');
    
cd ..


