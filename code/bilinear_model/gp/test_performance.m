function performance = test_performance(matname, normalization, plotting)
    if nargin < 2
        normalization = true;
    end
    if nargin < 3
        plotting = false;
    end
    
    GPMODEL = load(matname);
    
    % Load test data
    [xTest, yTest] = gp_load_data('test', GPMODEL.dist_type);
    signals = gp_create_signals(xTest, yTest, normalization);
    
    signals_model = SignalsModel(signals);
    signals_model.setInputs(GPMODEL.input_model);
    signals_model.setOutput('temp');
    
    % Obtain all IO data
    [X_test, Y_test] = signals_model.getIOVectors([], 'removeNaNs', true);
    
    [Ytest_mean, Ytest_var] = gp(GPMODEL.hyp, @infExact, GPMODEL.meanfunc, GPMODEL.covfunc, GPMODEL.likfunc, GPMODEL.X, GPMODEL.Y, X_test);
    
    % Denormalize output
    output_norm = signals.getSignalNormalization('temp');
    if ~isempty(output_norm)
        Y_test = SignalsValues.postNorm(Y_test, output_norm(1), output_norm(2));
        Ytest_mean = SignalsValues.postNorm(Ytest_mean, output_norm(1), output_norm(2));
        Ytest_var = SignalsValues.postNormVar(Ytest_var, output_norm(1), output_norm(2));
    end
    
    [performance.ae, performance.se, performance.lpd, performance.mrse, performance.smse, performance.msll] = ...
        loss(Y_test, Ytest_mean, Ytest_var);
    
    res = Ytest_mean - Y_test;
    ymean = mean(Y_test);
    performance.r2 = 1 - sum(res.^2)/sum((Y_test-ymean).^2);
    
    save([matname, '-test'], 'performance', 'Ytest_mean', 'Ytest_var', 'Y_test');
    
    if plotting
        h = figure;
        plotgp(h, (0:length(Y_test)-1)', Y_test, Ytest_mean, sqrt(Ytest_var));
        title(sprintf('%s - %s - %d samples', GPMODEL.covfuncshort, GPMODEL.meanfuncshort, GPMODEL.Ntrain));
    end
    
end