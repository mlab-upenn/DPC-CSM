
%% forest

plot_settings;
figure; hold on; grid on; box on;
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2.5];

startDate = datenum('05-20-2007');
endDate = datenum('05-26-2007');
t = linspace(startDate,endDate,168);

load('../results/dpcrf-validation.mat')
N = ctrlHzn;
h1 = plot(t, yPred(:,1), 'LineWidth', 2);
h3 = plot(t(N:end), yPred(1:end-N+1,N), 'LineWidth', 2);
h4 = plot(t, yTrue(:,1), 'LineWidth', 2);
hleg = legend([h1, h3, h4], 'k+1', 'k+N', 'ground truth');
set(hleg, 'location', 'south', 'Orientation', 'horizontal', 'box', 'on');
ylim([26 32]);

ylabel('room temp. $[\mathrm{^oC}]$');
xlabel('mm/dd [-]');
datetick('x','mm/dd')

%% tree

plot_settings;
figure; hold on; grid on; box on;
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2.5];

startDate = datenum('05-20-2007');
endDate = datenum('05-26-2007');
t = linspace(startDate,endDate,168);

load('../results/dpcrf-validation.mat')
N = ctrlHzn;
h1 = plot(t, yPred(:,1), 'LineWidth', 2);
h3 = plot(t(N:end), yPred(1:end-N+1,N), 'LineWidth', 2);
h4 = plot(t, yTrue(:,1), 'LineWidth', 2);
hleg = legend([h1, h3, h4], 'k+1', 'k+N', 'ground truth');
set(hleg, 'location', 'south', 'Orientation', 'horizontal', 'box', 'on');
ylim([26 32]);

ylabel('room temp. $[\mathrm{^oC}]$');
xlabel('mm/dd [-]');
datetick('x','mm/dd')
