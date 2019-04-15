
%% first step prediction

plot_settings;
figure; hold on; %grid on; box on;
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2];

startDate = datenum('05-20-2007');
endDate = datenum('05-27-2007');
t = linspace(startDate,endDate,168);

load('../results/dpcrf-validation.mat')
h4 = plot(t, yTrue(:,1), '-', 'LineWidth', 1);
h1 = plot(t, yPred(:,1), '-','LineWidth', 1);
yPredrf = yPred(:,1);
load('../results/dpcrt-validation.mat')
h3 = plot(t, yPred(:,1), '-', 'LineWidth', 0.5);
yTrue = yTrue(:,1);
yPredrt = yPred(:,1);
hleg = legend([h4, h1, h3], 'data', 'forest', 'tree');
set(hleg, 'location', 'southeast', 'Orientation', 'horizontal', 'box', 'off');
ylim([26 32]);

csvwrite('../results/validation1.csv', [yTrue, yPredrt, yPredrf]);

ylabel('$\mathsf{Y}_{\mathrm{k+1|k}}$ $[\mathrm{^oC}]$');
% xlabel('mm/dd');
datetick('x','mm/dd')

print('-depsc2', '-r600', '../latex/validation-s1.eps')

%% 6 hour ahead prediction

plot_settings;
figure; hold on; %grid on; box on;
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2];

load('../results/dpcrf-validation.mat')
h4 = plot(t, yTrue(:,ctrlHzn), '-', 'LineWidth', 1);
h1 = plot(t, yPred(:,ctrlHzn), '-', 'LineWidth', 1);
yPredrf = yPred(:,ctrlHzn);
load('../results/dpcrt-validation.mat')
h3 = plot(t, yPred(:,ctrlHzn), '-', 'LineWidth', 0.5);
yTrue = yTrue(:,ctrlHzn);
yPredrt = yPred(:,ctrlHzn);
% hleg = legend([h4, h1, h3], 'ground truth', 'forest', 'tree');
% set(hleg, 'location', 'southeast', 'Orientation', 'horizontal', 'box', 'on');
ylim([26 32]);

csvwrite('../results/validation6.csv', [yTrue, yPredrt, yPredrf]);

ylabel('$\mathsf{Y}_{\mathrm{k+6|k}}$ $[\mathrm{^oC}]$');
xlabel('mm/dd');
datetick('x','mm/dd')

print('-depsc2', '-r600', '../latex/validation-s6.eps')