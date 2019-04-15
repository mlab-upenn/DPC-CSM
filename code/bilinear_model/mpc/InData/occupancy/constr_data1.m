%% constr_data.m
% construct data out of actelion measurements
% divided into two parts:
% part 
%   a) get occupancy right
%   b) get electrical energy use correct
% v1: small modifications, with std deviation

clc; clear all; 
% close all
load('actelion_measurements_aprToAug.mat')

%% a) occupancy realization in 15min intervals
% collect occupancy data, all data of high quality
% rooms [B3 B4 B5 B7 B8 B9 B10 B11]
RoomPresence = nan(length(ctrllog.meas),8);

for ii = 1 : length(ctrllog.meas)
    RoomPresence(ii,:) = [  ctrllog.meas{ii}.data(5)
                            ctrllog.meas{ii}.data(8)
                            ctrllog.meas{ii}.data(11)
                            ctrllog.meas{ii}.data(17)
                            ctrllog.meas{ii}.data(20)
                            ctrllog.meas{ii}.data(23)
                            ctrllog.meas{ii}.data(26)
                            ctrllog.meas{ii}.data(29)   ]';
end

% figure(); hold on;
% title('total number of occupants on 2nd floor over single rooms')
% plot(RoomPresence_tot,'b','LineWidth',2)
% legend('tot occupants')
% hold off

%% do for individual rooms
start = 830;    % get start correct for 14 weeks, good start: 830
RoomPresence_sorted_ind = RoomPresence(start:end,:);    % individual rooms
RoomPresence_sorted_ind(end-10:end,:) = []; % correct weird data at end
toExtend = 13*7*24*4 - length(RoomPresence_sorted_ind); % extend to full 13 weeks
RoomPresence_sorted_ind = [RoomPresence_sorted_ind ; zeros(toExtend,8)];

% figure(); hold on;
% plot(RoomPresence_sorted_ind(:,2),'b','LineWidth',3)
% hold off

% change from 15min period to 1hr period
RoomPresence_sorted_ind_hrs = nan(13*7*24,8);
for ii = 1 : 13*7*24
    RoomPresence_sorted_ind_hrs(ii,:) = mean(RoomPresence_sorted_ind(1+(ii-1)*4:ii*4,:),1);
end

for ii = 1 : 8
    figure(); hold on;
    plot(RoomPresence_sorted_ind_hrs(:,ii),'b','Linewidth',2);
    title(['Hourly Presence Profile of room ' num2str(ii)])
    hold off
end


%% do for average of all rooms
RoomPresence_sorted = mean(RoomPresence_sorted_ind,2);   % average over all rooms
RoomPresence_sortedStd = std(RoomPresence_sorted_ind,0,2);

figure(); hold on;
stairs(RoomPresence_sorted,'b','LineWidth',2)
stairs(RoomPresence_sortedStd,':r','LineWidth',2)
legend('mean over 8 rooms','std deviation')
hold off;

% change from 15min to 1hr samples
RoomPresence_sorted_hrs = nan(13*7*24,1);
for ii = 1 : 13*7*24
    RoomPresence_sorted_hrs(ii) = mean(RoomPresence_sorted(1+(ii-1)*4 : ii*4));
end

% figure(); hold on;
% title('average occupancy profile 1 hr sample')
% stairs(RoomPresence_sorted_hrs,'b','LineWidth',2)
% hold off

% divide RoomPresence into weeks
% RoomPresence_sorted_hrs_div = time x 13
RoomPresence_sorted_hrs_div = nan(7*24,13);
for ii = 1 : 13
    RoomPresence_sorted_hrs_div(:,ii) = RoomPresence_sorted_hrs(1+(ii-1)*7*24:ii*7*24);
end

% pack eight rooms and 13 weeks together
RoomPresence_sorted_hrs_div_tot = mean(RoomPresence_sorted_hrs_div,2);

%% take the average of the week
% should be used as occupancy prediction as internal gains
RoomPresence_sorted_hrs_div_sched = RoomPresence_sorted_hrs_div_tot;
RoomPresence_sorted_hrs_div_constr = round(RoomPresence_sorted_hrs_div_sched*2);
SIA = [0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0]';
SIA = repmat(SIA,[5,1]);
SIA = [SIA ; zeros(48,1)];

% schedule needs a bit postprocessing for Friday
RoomPresence_sorted_hrs_div_constr(107) = 1;

figure(); hold on;
title('Average number of People per Room per Week');
stairs(RoomPresence_sorted_hrs_div_constr,'b','LineWidth',2)
stairs(SIA,'--c','LineWidth',2)
stairs(RoomPresence_sorted_hrs_div_sched,':r','LineWidth',2)
legend('Constraints for Actelion', 'Constraints SIA', 'Schedule: Avg number of people')
hold off

%% collect the working days together (ignore) weekends and take mean
RoomPresence_1day = nan(24,5);
for ii = 1 : 5
    RoomPresence_1day(:,ii) = RoomPresence_sorted_hrs_div_sched(1+(ii-1)*24:ii*24);
end

RoomPresence_1dayStd = std(RoomPresence_1day,0,2);  % compute std
RoomPresence_1day = mean(RoomPresence_1day,2);

% assume when there will be people
RoomPresence_on = RoomPresence_1day > 0.1;
RoomPresence_avg = sum(RoomPresence_1day(RoomPresence_1day > 0.2))/sum(RoomPresence_1day > 0.2);
RoomPresence_avg = RoomPresence_avg * RoomPresence_on;
RoomPresence_avg = 0.3*RoomPresence_on;

figure(); hold on;
title('Empirical Occupancy of week')
stairs(RoomPresence_1day,'b','LineWidth',2)
stairs(RoomPresence_1dayStd,':r','LineWidth',2)
stairs(RoomPresence_on,'g','LineWidth',2)
stairs(RoomPresence_avg,'-c','LineWidth',2)
legend('average room presence','std room presence','for setback','for prediction')
hold off

%% assuming two cases, see how distribution vary w.r.t 
% 1. RoomPresence_1day
% 2. RoomPresence_avg
% on a daily basis



%% save prediction for in 0-1 form
ig_act_constr_1we = [repmat(RoomPresence_on,[5,1]) ; zeros(2*24,1)];
% save('ig_act_constr_1we.mat','ig_act_constr_1we');

figure(); hold on;
title('for constraint setbacks')
plot(ig_act_constr_1we,'LineWidth',2)
hold off

% prediction
ig_act_pred_1we = RoomPresence_sorted_hrs_div_sched;
% save('ig_act_pred_1we.mat','ig_act_pred_1we')

figure(); hold on;
title('for prediction')
plot(ig_act_pred_1we,'LineWidth',2)
hold off


%% build realizations out of eight columns of 'RoomPresence_sorted_ind_hrs'
% take cols 1, 2, 3, 8
ig_act_real = [RoomPresence_sorted_ind_hrs(:,1)
               RoomPresence_sorted_ind_hrs(:,2)
               RoomPresence_sorted_ind_hrs(:,3)
               RoomPresence_sorted_ind_hrs(:,4)];
% save('ig_act_real.mat','ig_act_real')









