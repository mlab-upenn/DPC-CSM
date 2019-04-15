clear
clc

%% model training

% set params
orderAR = 6;
ctrlHzn = 6;
dist_type = 'est';

[xTrain_d, xTrain_c, yTrain] = load_data('train', orderAR, ctrlHzn, dist_type);

leafsize = 5;
catcol = 8*ctrlHzn+orderAR+(1:2);

% separate regression tree for each output
models = cell(1,ctrlHzn);
for idm = 1:ctrlHzn
    model = fitrtree(xTrain_d, yTrain(:,idm), 'categoricalpredictors', catcol,...
        'MinLeafSize', leafsize+4*idm);
    models{idm} = model;
    save(['../results/rt-' dist_type '-step' num2str(idm) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat'], 'model');
    disp(idm);
end

%% fit linear models in leaves

leafmodels = cell(1,ctrlHzn);
for idm = 1:ctrlHzn
    linmodel = train_linearmodel_in_leaves(models{idm}, xTrain_c, yTrain(:,idm), idm);
    leafmodels{idm} = linmodel;
    save(['../results/lm-rt-' dist_type '-step' num2str(idm) '-ctrlHzn' num2str(ctrlHzn) '-orderAR' num2str(orderAR) '.mat'], 'linmodel', '-v7.3');
end
