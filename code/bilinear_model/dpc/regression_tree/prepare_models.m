
models = cell(1,ctrlHzn);
leafmodels = cell(1,ctrlHzn);

disp('loading trees and leaf models...');

for idm = 1:ctrlHzn
    
    disp(['tree ' num2str(idm) ' of ' num2str(ctrlHzn)]);
    load(['results/rt-' data_type '-step' num2str(idm) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat'])
    models{idm} = model;
    load(['results/lm-rt-' data_type '-step' num2str(idm) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat'])
    leafmodels{idm} = linmodel;
    
end