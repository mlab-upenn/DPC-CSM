function [features, outputs] = gp_load_data(argStr, estType)
% Load data for GP models.
    
    if strcmp(argStr, 'train')
        if strcmp(estType, 'orig')
            DATA = load('../data/train-orig-filtered.mat');
        elseif strcmp(estType, 'est')
            DATA = load('../data/train-est-filtered.mat');    
        end
    elseif strcmp(argStr, 'test')
        if strcmp(estType, 'orig')
            DATA = load('../data/test-orig-filtered.mat');
        elseif strcmp(estType, 'est')
            DATA = load('../data/test-est-filtered.mat');
        end
    else
        error('input string should be either ''train'' or ''test'' ');
    end

    % Column 3 of disturbance inputs is constant 0, so remove it
    DATA.d(3,:) = [];

    % Create the fields in the resulting structures to hold different time series
    N = size(DATA.d, 1);
    features = struct;
    for k = 1:N
        features.(sprintf('d%d', k)) = DATA.d(k,:);
    end
    features.tod = DATA.proxy(1,:);
    features.dow = DATA.proxy(2,:);

    % Control inputs
    features.blind = DATA.u(1,:);
    features.light = DATA.u(2,:);
    features.heating = DATA.u(3,:);
    features.cooling = DATA.u(4,:);

    % Output = room temperature
    outputs = DATA.y(1,:);
end