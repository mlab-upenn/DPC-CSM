
models = cell(1,ctrlHzn);
leafmodels = cell(1,ctrlHzn);

disp('loading forests and leaf models...');

for idm = 1:ctrlHzn
    
    disp(['forest ' num2str(idm) ' of ' num2str(ctrlHzn)]);
    load(['results/rf-' data_type '-step' num2str(idm) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat'])
    models{idm} = model;
    load(['results/lm-rf-' data_type '-step' num2str(idm) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat'])
    leafmodels{idm} = linmodel;
    
end
