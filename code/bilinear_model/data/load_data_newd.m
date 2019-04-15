% idm is the tree number on the horizon we want to load disturbance for

function [features_d, features_c, outputs] = load_data_newd(argStr, orderAR, ctrlHzn, estType, idm)

if strcmp(argStr, 'train')
    if strcmp(estType, 'orig')
        load('train-orig-filtered.mat');
    elseif strcmp(estType, 'est')
        load('train-est-filtered.mat');    
    end
elseif strcmp(argStr, 'test')
    if strcmp(estType, 'orig')
        load('test-orig-filtered.mat');
    elseif strcmp(estType, 'est')
        load('test-est-filtered.mat');
    end
else
    error('input string should be either ''train'' or ''test'' ');
end

features_d = [lagmatrix(d', 0-idm+1:1:ctrlHzn-1-idm+1), lagmatrix(y(1,:)', 0:1:orderAR-1), proxy']; %#ok<*NODEF>
features_c = lagmatrix(u', -(0:1:ctrlHzn-1));
outputs = lagmatrix(y(1,:)', -(1:1:ctrlHzn));

features_d(1:orderAR-1,:) = [];
features_c(1:orderAR-1,:) = [];
outputs(1:orderAR-1,:) = [];
features_d(end-ctrlHzn:end,:) = [];
features_c(end-ctrlHzn:end,:) = [];
outputs(end-ctrlHzn:end,:) = [];