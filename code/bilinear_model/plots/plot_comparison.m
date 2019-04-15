plot_settings;

close all;

data_type = 'est';      % est/orig
dpc_type = 'tree';    % tree/forest
orderAR = 6;
ctrlHzn = 6;

tstart = 744-3*24+1;
tstop = 744+3*24;  
steps = tstart:tstop;
nsteps = length(steps);
DRbegin = 21;
DRend = 44;
ecost = 1;

startDate = datenum('01-29-2007');
endDate = datenum('02-4-2007');
t = linspace(startDate,endDate,144);
may1 = datenum('05-01-2007');
may2 = datenum('05-02-2007');
may3 = datenum('05-03-2007');

fig1 = figure(1); hold on;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 6 2.2];

fig2 = figure(2); hold on;
fig2.PaperUnits = 'inches';
fig2.PaperPosition = [0 0 6 2];

fig3 = figure(3); hold on;
fig3.PaperUnits = 'inches';
fig3.PaperPosition = [0 0 6 1.6];

fig4 = figure(4); hold on;
fig4.PaperUnits = 'inches';
fig4.PaperPosition = [0 0 6 1.3];

fig5 = figure(5); hold on;
fig5.PaperUnits = 'inches';
fig5.PaperPosition = [0 0 6 2];

% mpc
load(['../results/mpc-' data_type '-cost' num2str(ecost) '-start' num2str(tstart) '-stop' num2str(tstop) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat']);

figure(1);
h11 = plot(t, yvec(1,steps), 'LineWidth', 1);
% ymax_all(ymax_all(:,1)>35,:) = 27;
% ymin_all(ymin_all(:,1)<10,:) = 19;
plot(t, ymax_all(steps,1), '--k', 'LineWidth', 1)
plot(t, ymin_all(steps,1), '--k', 'LineWidth', 1)
ylim([17,28])
ylabel('room temp. $[\mathrm{^oC}]$');
xlabel('mm/dd');
datetick('x','mm/dd');
set(gca,'XTickLabel',{'01/29','01/30','01/31','05/01','05/02','05/03','05/04'})

figure(2);
h21 = plot(t,uvec(3,steps), 'LineWidth', 1);
plot(t,umax_all(steps,3), '--k', 'LineWidth', 1)
plot(t,umin_all(steps,3), '--k', 'LineWidth', 1)
ylabel('cooling factor $[-]$');
ylim([0,1.5])
% xlabel('mm/dd');
datetick('x','mm/dd')
set(gca,'XTickLabel',{'01/29','01/30','01/31','05/01','05/02','05/03','05/04'})

figure(3);
h31 = plot(t,uvec(4,steps), 'LineWidth', 1);
plot(t,umax_all(steps,4), '--k', 'LineWidth', 1);
plot(t,umin_all(steps,4), '--k', 'LineWidth', 1);
ylabel('heat $[\mathrm{W/m^2}]$');
ylim([0,24])
% xlabel('mm/dd');
datetick('x','mm/dd')
set(gca,'XTickLabel',{'01/29','01/30','01/31','05/01','05/02','05/03','05/04'})

figure(4);
h41 = plot(t,cumsum(fvalvec(1,steps)), 'LineWidth', 1);
ylabel('cum. cost');
ylim([0,6000])
% xlabel('mm/dd');
datetick('x','mm/dd')
set(gca,'XTickLabel',{'01/29','01/30','01/31','05/01','05/02','05/03','05/04'})

% dpcrf
load(['../results/dpcrf-' data_type '-cost' num2str(ecost) '-start' num2str(tstart) '-stop' num2str(tstop) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat']);
figure(1);
h12 = plot(t,yvec(1,steps), 'LineWidth', 1);
figure(2);
h22 = plot(t,uvec(3,steps), 'LineWidth', 1);
figure(3);
h32 = plot(t,uvec(4,steps), 'LineWidth', 1);
figure(4);
h42 = plot(t,cumsum(fvalvec(1,steps)), 'LineWidth', 1);

% dpcrt
load(['../results/dpcrt-' data_type '-cost' num2str(ecost) '-start' num2str(tstart) '-stop' num2str(tstop) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat']);
figure(1);
h13 = plot(t,yvec(1,steps), '-', 'LineWidth', 0.5);
% hleg = legend([h11, h12, h13], 'DPC-RT', 'DPC-En', 'MPC');
% set(hleg, 'location', 'north', 'Orientation', 'horizontal', 'box', 'on');
print('-depsc2', '-r600', '../latex/state.eps')

figure(2);
h23 = plot(t,uvec(3,steps), '-', 'LineWidth', 0.5);
hleg = legend([h21, h22, h23], 'MPC', 'DPC-En', 'DPC-RT');
set(hleg, 'location', 'north', 'Orientation', 'horizontal', 'box', 'off');
print('-depsc2', '-r600', '../latex/input3.eps')

figure(3);
h33 = plot(t,uvec(4,steps), '-', 'LineWidth', 0.5);
% hleg = legend([h31, h32, h33], 'DPC-RT', 'DPC-En', 'MPC');
% set(hleg, 'location', 'north', 'Orientation', 'horizontal', 'box', 'on');
print('-depsc2', '-r600', '../latex/input4.eps')

figure(4);
h43 = plot(t,cumsum(fvalvec(1,steps)), 'LineWidth', 1);
print('-depsc2', '-r600', '../latex/cumsumcost.eps')

figure(5);
h1 = plot(t,dvec(2,steps), 'LineWidth', 1);
h2 = plot(t,dvec(6,steps), 'LineWidth', 1);
h3 = plot(t,dvec(7,steps), 'LineWidth', 1);
ylim([-5,60]);
% xlabel('mm/dd');
datetick('x','mm/dd')
set(gca,'XTickLabel',{'01/29','01/30','01/31','05/01','05/02','05/03','05/04'})
hleg = legend([h1,h2,h3], 'solar gain \ $[\mathrm{W/m^2}]$', 'equip. gain $[\mathrm{W/m^2}]$', 'dry-bulb temp. $[\mathrm{^oC}]$');
text(t(10),53,'\textbf{see colors in digital version}');
set(hleg, 'location', 'northeast', 'Orientation', 'vertical', 'box', 'off');
print('-depsc2', '-r600', '../latex/disturbances.eps')