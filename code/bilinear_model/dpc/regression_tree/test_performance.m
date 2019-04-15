orderAR = 6;
ctrlHzn = 6;
data_type = 'est';

prepare_models;

[xTest_d, xTest_c, yTest] = load_data('test', orderAR, ctrlHzn, data_type);

n = 744+19*24+1:744+26*24;
xTest_d = xTest_d(n,:);
xTest_c = xTest_c(n,:);
yTest = yTest(n,:);

% calculate root mean square error with our definition
yTrue = yTest;
yPred = zeros(size(xTest_d,1),ctrlHzn);

% find linear models
for idn = 1:size(xTest_d,1)
    for idm = 1:ctrlHzn
        coeff = find_linearmodel_in_leaves(models{idm}, leafmodels{idm}, xTest_d(idn,:));
        yPred(idn,idm) = [1, xTest_c(idn,1:4*idm)]*coeff;
    end
    disp(idn);
end

nrmseTest = sqrt(mean((yPred(ctrlHzn+1:end,:)-yTrue(ctrlHzn+1:end,:)).^2,1))./mean(yTrue(ctrlHzn+1:end,:),1);
fprintf('nrmse on test data using our definition = %f\n', min(nrmseTest));
% plot_ml_results(yTrue, yPred', 'true-predicted', ctrlHzn);
plot_ml_results(yTrue, yPred', 'predictive', ctrlHzn)

save('../../results/dpcrt-validation.mat','n','xTest_d','xTest_c','yTest','yPred','yTrue','orderAR','ctrlHzn','data_type');