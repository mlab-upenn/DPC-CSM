%% constr_data.m
% construct data out of actelion measurements
% divided into two parts:
% part 
%   a) get occupancy right
%   b) get electrical energy use correct
% v1: small modifications, with std deviation
% v2: use Rooms 3, 6, 7, 8 to compute std, approx sched, and dist
%     use other Rooms as realizations of one year

clc; clear all; 
% close all
load('actelion_measurements_aprToAug.mat')

% for ii = 1 : 25
    tmp = randperm(8);
    useRooms = tmp(1:4);
    useRooms = [3 6 7 8];

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

    %% do for individual rooms
    start = 830;    % get start correct for 14 weeks, good start: 830
    RoomPresence_sorted_ind = RoomPresence(start:end,:);    % individual rooms
    RoomPresence_sorted_ind(end-10:end,:) = []; % correct weird data at end
    toExtend = 13*7*24*4 - length(RoomPresence_sorted_ind); % extend to full 13 weeks
    RoomPresence_sorted_ind = [RoomPresence_sorted_ind ; zeros(toExtend,8)];

    % change from 15min period to 1hr period
    RoomPresence_sorted_ind_hrs = nan(13*7*24,8);
    for ii = 1 : 13*7*24
        RoomPresence_sorted_ind_hrs(ii,:) = mean(RoomPresence_sorted_ind(1+(ii-1)*4:ii*4,:),1);
    end

    % divide rooms into two groups
    RoomPresence_pred = [];     % used for predictions
    RoomPresence_real = [];     % used as realizations in simulations

    for ii = 1 : 8
        if sum(useRooms == ii)
            RoomPresence_pred = [RoomPresence_pred RoomPresence_sorted_ind_hrs(:,ii)];
        else
            RoomPresence_real = [RoomPresence_real RoomPresence_sorted_ind_hrs(:,ii)];
        end
    end

    % divide into pieces of weeks
    RoomPresence_pred_weeks = [];
    RoomPresence_real_weeks = [];
    for ii = 1 : 13
        RoomPresence_pred_weeks = [RoomPresence_pred_weeks RoomPresence_pred(1+(ii-1)*7*24:ii*7*24,:)];
        RoomPresence_real_weeks = [RoomPresence_real_weeks RoomPresence_real(1+(ii-1)*7*24:ii*7*24,:)];
    end

    RoomPresence_pred_weeks_mean = mean(RoomPresence_pred_weeks,2);
    RoomPresence_real_weeks_mean = mean(RoomPresence_real_weeks,2);
    RoomPresence_pred_weeks_std = std(RoomPresence_pred_weeks,0,2);
    RoomPresence_real_weeks_std = std(RoomPresence_real_weeks,0,2);

    figure(); hold on;
    title(['Average Occupancy ' num2str(useRooms)])
    plot(RoomPresence_pred_weeks_mean,'b','LineWidth',2)
    plot(RoomPresence_real_weeks_mean,'--r','LineWidth',2)
    legend('mean prediction','mean realization')

    figure(); hold on;
    title(['Std Occupancy ' num2str(useRooms)])
    plot(RoomPresence_pred_weeks_std,'c','LineWidth',2)
    plot(RoomPresence_real_weeks_std,'--m','LineWidth',2)
    legend('std prediction','std realization')
    
% end

%% prepare to save prediction
ig_act_pred_1we = RoomPresence_pred_weeks_mean;
% save('ig_act_pred_1we.mat','ig_act_pred_1we')

% save realization
ig_act_real = reshape(RoomPresence_real,[2184*4,1]);
ig_act_real = [ig_act_real ; ig_act_real(1:48)];

% save('ig_act_real.mat','ig_act_real')

%% overview of what happens
igPred = repmat(ig_act_pred_1we,[52,1]);
load('ig_act_constr_1we.mat');
igConstr = repmat(ig_act_constr_1we,[52,1]);

figure(); hold on;
title('occupancy data')
plot(igConstr,'k')
plot(igPred,'b')
plot(ig_act_real,'r')
plot(RoomPresence_pred_weeks_std,'m')
legend('comfort constraints','Predicted IG','Realization IG','Std pred weeks')
hold off


%% using past predictions, get std in distribution during 'on'-time
std_avg = RoomPresence_pred_weeks_std .* ig_act_constr_1we;
std_avg = std_avg(std_avg >0);
std_avg = mean(std_avg);
disp(['average std dev: ' num2str(std_avg)]);

%% plot all curves
figure(); hold on;
title('RoomPresence of 4 rooms over 13 weeks')
pert = 0.01*(rand(size(RoomPresence_pred_weeks))-0.5);
RoomPresence_pred_weeks_pert = RoomPresence_pred_weeks + pert;
for ii = 1 : size(RoomPresence_pred_weeks,2)
    plot(RoomPresence_pred_weeks_pert(:,ii),'o', 'Color', [0.7,0.7,0.7]);
end
hold off


%% construct samples of deviations from old data
ig_act_realPred = reshape(RoomPresence_pred,[2184*4,1]);    % use old Data
ig_act_pred_1year = repmat(ig_act_pred_1we,[52,1]);
errors_1year = ig_act_realPred - ig_act_pred_1year;

% sort errors according to 
%  1. Day: 1,...,7
%  2. time: 0,..,23
% this ways, we can avoid pred+dist that are >1 or <0

% first divide errors into hours
errDeltaIG_h = cell(24,1);
for jj = 1 : 24
    errDeltaIG_h{jj} = [];
    for ii = jj : 24 : length(errors_1year)
        errDeltaIG_h{jj} = [errDeltaIG_h{jj}
                            errors_1year(ii)];
    end
end

% divide things into 7 days
errDeltaIG_h_d = cell(24,7);
for jj = 1 : 24     % iterate over 24 hours
    for k = 1 : 7   % iterate over 7 days
        errDeltaIG_h_d{jj,k} = errDeltaIG_h{jj}(k:7:end);
    end
end

save('errDeltaIG_h_d.mat','errDeltaIG_h_d');


%% throw all errors together into one bin
% extract errors only during daytime
extract = repmat(ig_act_constr_1we,[52,1]);
errDeltaIG_1year = errors_1year(extract > 0);

save('errDeltaIG_1year.mat','errDeltaIG_1year');

[n x] = hist(errDeltaIG_1year,100);
n = n/sum(n);
figure(); hold on;
title('Histogram of prediction error of 1 year, daytime')
plot(x,n)
hold off

%% plot for report
figure(); hold on;
title('Internal Gain Persons one Week', 'FontSize', 18)
stairs(ig_act_real(1:7*24), 'b', 'LineWidth',2)
stairs(ig_act_pred_1we,'--r','LineWidth',2)
h_legend = legend('Realization', 'Average / Prediction')
xlabel('time steps [h]', 'FontSize', 18)
ylabel('Persons', 'FontSize', 18)
set(h_legend, 'FontSize', 16)
hold off

