%%  actelion_measurements
%   extract and process data
%   Date: 05. August
%   Author: Xiaojing Zhang

% assumptions: Electrical Energy consumption as follows per person:
%   a) computer + display: 280 + 20 = 300W
%   b) electrical energy used (examples):
%        i.     8.2 sensors on day: 607 Watt
%        ii.    8.65 sensors on day: 455 Watt
%        iii.   5 sensors on day: 480 Watt
% insights: 
%   1. for energy consumption, use measured total energy minus light energy
%   2. take working days as 5 days + 2 days, not the strange days

clear all; close all; clc;

load('actelion_measurements_aprToAug.mat')

%% process data

figure(); hold on;
title('time')
plot(ctrllog.dtnum)
hold off


%% collect global radiation (74)th element
% couple of errors

RGS = nan(length(ctrllog.meas),1);

for ii = 1 : length(ctrllog.meas)
    RGS(ii) = ctrllog.meas{ii}.data(74);
    if ~ctrllog.meas{ii}.quality(74)
        disp(['quality error in RGS at time ' num2str(ii)])
    end
end


figure(); hold on;
title('RGS')
plot(RGS)
% plot(ctrllog.dtnum(1:length(ctrllog.meas)),RGS,'LineWidth',2)
% datetick('x',31,'keepticks')
legend('RGS')
hold off

%% collect ElEnergy counters
% no errors
% Kanal [S E NO N NW W total] in cumulative total energy 

ElEnergy = nan(length(ctrllog.meas),7);     % kWh

for ii = 1 : length(ctrllog.meas)
    ElEnergy(ii,:) = ctrllog.meas{ii}.data(59:65)';
    if sum(ctrllog.meas{ii}.quality(59:65)) ~= 7
        disp(['quality error in ElEnergy at time ' num2str(ii)])
    end
end

ElEnergy_tot = sum(ElEnergy(:,1:6),2);

figure(); hold on;
title('total commutative Electrical Energy [kWh]')
% plot(ctrllog.dtnum(1:length(ctrllog.meas)), ElEnergy(:,end),'b','LineWidth',2);
% plot(ctrllog.dtnum(1:length(ctrllog.meas)), ElEnergy_tot,'--r','LineWidth',2);
plot(ElEnergy(:,end),'b','LineWidth',2);
plot(ElEnergy_tot,'--r','LineWidth',2);
legend('automatically from counter','manually added');
% datetick('x',31,'keepticks')
hold off;

%% collect ElPower 2OG
% no error
% Kanal [S E NO N NW W total] in cumulative total energy 

ElPower = nan(length(ctrllog.meas),7);     % W

for ii = 1 : length(ctrllog.meas)
    ElPower(ii,:) = ctrllog.meas{ii}.data(103:109)';
    if sum(ctrllog.meas{ii}.quality(103:109)) ~= 7
        disp(['quality error in ElPower at time ' num2str(ii)])
    end
end

ElPower_tot = sum(ElPower(:,1:6),2);

figure(); hold on;
title('total Electrical Power [W]');
% plot(ElPower(:,end),'b','LineWidth',2);
% plot(ElPower_tot,'--r','LineWidth',2);
plot(ElPower(:,end),'b','LineWidth',2);
plot(ElPower_tot,'--r','LineWidth',2);
legend('autom from counter', 'manually added')
hold off

%% collect Room Presence together
% all high quality
% Kanal Rooms [B1 B3 B4 B5 B6 B7 B8 B9 B10 B11 B13 B14 B15 B17 B19 B21 B23 B24 B26]
RoomPresence = nan(length(ctrllog.meas),19);

for ii = 1 : length(ctrllog.meas)
    RoomPresence(ii,:) = ctrllog.meas{ii}.data(2:3:56)';
    if sum(ctrllog.meas{ii}.quality(2:3:56)) ~= 19
        disp(['quality error in RoomPresence at time ' num2str(ii)])
    end
end

RoomPresence_tot = sum(RoomPresence,2);
figure(); hold on;
title('total number of occupants on 2nd floor over 19 rooms')
% plot(ctrllog.dtnum(1:length(ctrllog.meas)),RoomPresence_tot,'b','LineWidth',2)
plot(RoomPresence_tot,'b','LineWidth',2)
legend('tot occupants')
% datetick('x','mmm','keepticks')
hold off

%% collect Electrical Energy consumption of Rooms, no data for room B26
% all high quality
% Kanal Rooms [B1 B3 B4 B5 B6 B7 B8 B9 B10 B11 B13 B14 B15 B17 B19 B21 B23 B24]
ElEnergyLight = nan(length(ctrllog.meas),18);   %[kWh]

for ii = 1 : length(ctrllog.meas)
    ElEnergyLight(ii,:) = ctrllog.meas{ii}.data(3:3:54)';
    if sum(ctrllog.meas{ii}.quality(3:3:54)) ~= 18
        disp(['quality error in RoomPresence at time ' num2str(ii)]);
    end
end

ElEnergyLight_tot = sum(ElEnergyLight,2);
figure(); hold on;
title('total added ElEnergyLight of 18 rooms')
plot(ElEnergyLight_tot,'b','LineWidth',2)
legend('tot ElEnergy Light consumption')
% datetick('x','mmm','keepticks')
hold off

%%
for roomNr = 1 : 19
    figure(); hold on;
    title(['number of occupants in Room ' num2str(roomNr)])
    plot(RoomPresence(:,roomNr))
    hold off
end

%% check quality
% not all rooms are used regularly
% quality drop in in solar radiation measurements

for ii=1 : length(ctrllog.meas)
    if sum(ctrllog.meas{ii}.quality) ~= length(ctrllog.meas{ii}.quality)
        disp(['error at time ' num2str(ii)])
    end
end

%% compute consumed power only between time 302 and 9200
ElEnergy_used_meas = ElEnergy(9200,end) - ElEnergy(302,end);    % kWh
ElEnergy_added_meas = ElEnergy_tot(9200) - ElEnergy_tot(302);     % kWh
ElPower_used_meas = sum(ElPower(302:9200,end))/4/1000;          % kWh
ElPower_added_meas = sum(ElPower_tot(302:9200,end))/4/1000;

disp(['ElEnergy_used_meas: ' num2str(ElEnergy_used_meas) ' kWh']);
disp(['ElEnergy_added_meas: ' num2str(ElEnergy_added_meas) ' kWh']);
disp(['ElPower_used_meas: ' num2str(ElPower_used_meas) ' kWh']);
disp(['ElPower_added_meas: ' num2str(ElPower_added_meas) ' kWh']);









