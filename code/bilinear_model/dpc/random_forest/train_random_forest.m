clear
clc

%% model training

% set params
orderAR = 6;
ctrlHzn = 6;

[xTrain_d, xTrain_c, yTrain] = load_data('train', orderAR, ctrlHzn, 'orig');

leafsize = 5;
catcol = 8*ctrlHzn+orderAR+(1:2);

% separate regression tree for each output
models = cell(1,ctrlHzn);
for idm = 1:ctrlHzn
    model = TreeBagger(400, xTrain_d, yTrain(:,idm),...
    'Method', 'regression', 'OOBPred', 'On', 'OOBVarImp', 'on',...
    'CategoricalPredictors', catcol, 'MinLeaf', leafsize+4*idm);
    models{idm} = model;
    save(['../results/rf-step' num2str(idm) '-N' num2str(ctrlHzn) '.mat'], 'model');
    disp(idm);
end

% clearvars -except model ctrlHzn orderAR xTrain_d xTrain_c yTrain
% save(['../results/models_random_forest_' num2str(ctrlHzn) '.mat'],'-v7.3')

%% fit linear models in leaves

leafmodels = cell(1,ctrlHzn);
for idm = 1:ctrlHzn
    linmodel = train_linearmodel_in_leaves(models{idm}, xTrain_c, yTrain(:,idm), idm);
    leafmodels{idm} = linmodel;
    save(['../results/lm-step' num2str(idm) '-N' num2str(ctrlHzn) '.mat'], 'linmodel', '-v7.3');
end

% clearvars -except leafmodels
% save(['../results/lin_random_forest_' num2str(ctrlHzn) '.mat'],'-v7.3')
