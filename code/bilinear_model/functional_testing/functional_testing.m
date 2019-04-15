%% load model

% already in discrete form
model = load('../data/bm.mat');
model = model.bm;

Ts = 3600; % sampling time

%% simulate model with random inputs

time = Ts*(0:1:8592-1);
Nu = 2; % number of time steps for which the input is constant

rng(0);
umin = [0      0          0           0         ]';
umax = [1      7.2        1           2*11.488  ]';
uRand = rand(4,length(time));
u = zeros(4,length(time));
for k = 1:length(time)
    u(:,k) = umin + (umax-umin).*uRand(:,floor(k/Nu)+1);
end
% u(1,:) = randi([0,1],[1,length(time)]);
% u(3,:) = randi([0,1],[1,length(time)]);

% load disturbance
load('../data/dist.mat')
isEst = 1;
if isEst
    dist_type = 'est';
    d = distEst';
else
    dist_type = 'orig';
    d = dist';
end

% initialize
x = zeros(12,length(time));
y = zeros(3,length(time));
x(:,1) = [22 22 22 22 22 22 22 22 22 20 15 10]';

cd('../model/')
% simulate the system
for k = 1:length(time)-1
    [x(:,k+1), y(:,k+1)] = simulate_model(model, x(:,k), d(:,k), u(:,k));
end

proxy = [mod((0:size(x,2)-1),24)'+1, mod(floor((0:size(x,2)-1)/24),7)'+1]';

cd('../functional_testing/')

%% plots

ifPlot = 1;

if ifPlot
    
    % disturbance plot 1
    figure;
    subplot(4,1,1); grid on;
    plot(time, d(1,:), 'LineWidth',1.5);
    legend('solar gain closed blinds')
%     axis([0 time(end) -20 40]);
    
    subplot(4,1,2); grid on;
    plot(time, d(2,:), 'LineWidth',1.5);
%     axis([0 time(end) -100 600]);
    legend('solar gain open blinds')
    
    subplot(4,1,3); grid on;
    plot(time, d(3,:), 'LineWidth',1.5);
%     axis([0 time(end) -100 1100]);
    legend('daylight illum. closed blinds')
    xlabel('time');
    
    subplot(4,1,4); grid on;
    plot(time, d(4,:), 'LineWidth',1.5);
%     axis([0 time(end) -100 1100]);
    legend('daylight illum. open blinds')
    xlabel('time');
    
    % disturbance plot 2
    figure;
    subplot(4,1,1); grid on;
    plot(time, d(5,:), 'LineWidth',1.5);
    legend('internal gain persons')
%     axis([0 time(end) -20 40]);
    
    subplot(4,1,2); grid on;
    plot(time, d(6,:), 'LineWidth',1.5);
%     axis([0 time(end) -100 600]);
    legend('internal gain equipment')
    
    subplot(4,1,3); grid on;
    plot(time, d(7,:), 'LineWidth',1.5);
%     axis([0 time(end) -100 1100]);
    legend('outside temp')
    xlabel('time');
    
    subplot(4,1,4); grid on;
    plot(time, d(8,:), 'LineWidth',1.5);
%     axis([0 time(end) -100 1100]);
    legend('free cooling temp')
    xlabel('time');
    
    figure; hold on; grid on;
    title('outputs')
    h1 = plot(time,y(1,:));
    h2 = plot(time,y(2,:));
    h3 = plot(time,y(3,:));
    set(h1, 'LineWidth',1.5);
    set(h2, 'LineWidth',1.5);
    set(h3, 'LineWidth',1.5);
    xlabel('time')
    legend([h1, h2, h3], 'room temp', 'room illum', 'ceiling temp')
%     axis([0 time(end) -20 40]);
    
    figure; hold on; grid on;
    title('inputs')
    plot(time, u(1,:), 'LineWidth',1.5);
    plot(time, u(2,:), 'LineWidth',1.5);
    plot(time, u(3,:), 'LineWidth',1.5);
    plot(time, u(4,:), 'LineWidth',1.5);
    xlabel('time')
    ylabel('inputs')
    legend('blind position', 'lighting', 'free cooling', 'radiator heating')
    
end

%% save results

save(['../data/data-' dist_type '.mat'], 'x', 'u', 'd', 'y', 'proxy');
csvwrite(['../../python/data/train_' dist_type '_x.csv'],x');
csvwrite(['../../python/data/train_' dist_type '_u.csv'],u');
csvwrite(['../../python/data/train_' dist_type '_d.csv'],d');
csvwrite(['../../python/data/train_' dist_type '_y.csv'],y');
